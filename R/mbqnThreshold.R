#library(ggplot2)
#library(reshape2)

#' @name mbqnGetThreshold 
#' @title Calculates the rank invariance threshold from which on MBQN should be used instead of 'classical' QN
#' @description Calculates the rank invariance threshold from which on MBQN should be used instead of 'classical' QN
#' @param mtx Matrix with samples in columns and features in rows
#' @param meanMedian Offset function for the MBQN calculation
#' @param plot Boolean values if logistic regression curves that are used to calculate intersection point should be plotted
#' @return threshold value
#' @export
mbqnGetThreshold <- function(mtx, meanMedian="mean", plot=TRUE) {
  mbqn.mtx <- mbqn(mtx,FUN = meanMedian)
  qn.mtx <- mbqn(mtx,FUN = NULL)
  
  pvalues_mbqn <- getPvalue(mtx, mbqn.mtx)
  pvalues_qn <- getPvalue(mtx, qn.mtx)
  
  low_thr <- 0.0
  res  <- mbqnGetNRIfeatures(mtx, low_thr = low_thr)

  nri_df <- data.frame(res$nri)
  nri_freq <- c()
  for (row in seq_len(nrow(nri_df))){
    freq <- nri_df[row,]$Freq
    nri_freq <- c(nri_freq, freq)
  }
  
  # Contain RI, p value and statistic
  combined_mbqn <- data.frame(nri_freq, pvalues_mbqn)
  combined_qn <- data.frame(nri_freq, pvalues_qn)
  
  thr = 0.01
  mbqnGetIntersect(combined_qn, combined_mbqn, thr, plot)
}

#' @name mbqnGetIntersect 
#' @title Helper function for mbqnGetThreshold
#' @description Helper function for mbqnGetThreshold
#' @param combined_qn Data frame containing RI, p value and statistic calculated for QN
#' @param combined_mbqn Data frame containing RI, p value and statistic calculated for MBQN
#' @param threshold Significance threshold for p value of Pitman-Morgan variance test
#' @param plot Boolean values if logistic regression curves that are used to calculate intersection point should be plotted
#' @return threshold value
#' @export
mbqnGetIntersect <- function(combined_qn, combined_mbqn, threshold, plot=TRUE){
  set.seed(2)
  thr <- 0.01
  combined_qn$sign <- ifelse(combined_qn$pvalue < thr, 1, 0)
  combined_mbqn$sign <- ifelse(combined_mbqn$pvalue < thr, 1, 0)
  
  model_qn <- stats::glm(sign ~ nri_freq, data = combined_qn, family = "binomial")
  model_mbqn <- stats::glm(sign ~ nri_freq, data = combined_mbqn, family = "binomial")
  
  if (model_qn$rank == 1){
    print("QN: GLM model rank is 1.")
  }
  
  if (model_mbqn$rank == 1){
    print("MBQN: GLM model rank is 1.")
  }
  
  z_value_qn <- summary(model_qn)$coefficients[2,3]
  sign_value_qn <- summary(model_qn)$coefficients[2,4]
  sign_value_qn <- oneSidedTest(sign_value_qn, z_value_qn)
  
  print(paste("Sign right-tailed QN: ", round(sign_value_qn, 3)))
  
  z_value_mbqn <- summary(model_mbqn)$coefficients[2,3]
  sign_value_mbqn <- summary(model_mbqn)$coefficients[2,4]
  sign_value_mbqn <- oneSidedTest(sign_value_mbqn, z_value_mbqn)
  
  print(paste("Sign right-tailed MBQN: ", round(sign_value_mbqn, 3)))
  
  # extract fitted probabilities
  predProbs_qn = stats::predict(model_qn, type="response")
  predProbs_mbqn = stats::predict(model_mbqn, type="response")
  
  f1 <- stats::approxfun(predProbs_qn - predProbs_mbqn, combined_qn$nri_freq, rule=2)
  
  # truncate to two decimal digits
  intersect <- truncateDecimals(f1(0)/100, digits=2)
  print(paste('Intersect: ',intersect))
  
  if (sign_value_qn < 0.01) { #0.001
    if (intersect == min(combined_qn$nri_freq)){
      print('QN Problem: Use MBQN for whole dataset.')
    } else {
      print(paste('QN Problem: Use MBQN. NRI-threshold =',intersect,'suggested.'))
    }
  }
  
  if (plot){
    ggplot2::theme_set(ggplot2::theme_minimal())
    plot_combined <- ggplot2::ggplot() +
      ggplot2::geom_point(ggplot2::aes(x=combined_qn$nri_freq, combined_qn$sign), color="#CC6666", alpha=0.3) +
      ggplot2::geom_line(ggplot2::aes(x=combined_qn$nri_freq, y=predProbs_qn), color="#CC6666") +
      ggplot2::geom_point(ggplot2::aes(x=combined_mbqn$nri_freq, combined_mbqn$sign), color="#66CC99", alpha=0.3) +
      ggplot2::geom_line(ggplot2::aes(x=combined_mbqn$nri_freq, y=predProbs_mbqn), color="#66CC99") +
      ggplot2::ylim(0,1) +
      ggplot2::xlim(0,100) +
      ggplot2::labs(
        x = "Frequency of RI proteins [%]",
        y = "Probability of significant p-values",
        title = paste("Threshold:", threshold, ", QN: red (", round(sign_value_qn, 3),"), MBQN: green (", round(sign_value_mbqn, 3), ")"))
    print(plot_combined)
    #ggsave(paste0("logistRegr_nriVsProbab_QNandMBQN_thresh_", threshold, ".pdf"), plot = plot_combined, width=10, height=8)
  }
  intersect
}


#library(PairedData)
#' @name getPvalue 
#' @title Calculates Pitman-Morgan variance test on two matrices
#' @description Calculates Pitman-Morgan variance test on two matrices
#' @param mtx1 Matrix with samples in columns and features in rows
#' @param mtx2 Matrix with samples in columns and features in rows
#' @return Data frame with p values and statistics
#' @export
getPvalue <- function(mtx1, mtx2){
  protnames <- c()
  pvalues <- c()
  statistics <- c()
  options("scipen"=100, "digits"=4)
  
  if (is.null(row.names(mtx1))){
    rownames(mtx1) <- 1:nrow(mtx1)
  }
  
  for (i in seq_len(nrow(mtx1))){
    protnames <- c(protnames, row.names(mtx1)[i])
    pvalues <- c(pvalues, PairedData::Var.test(mtx1[i,], mtx2[i,], paired=TRUE)$p.value)
    statistics <- c(statistics, PairedData::Var.test(mtx1[i,], mtx2[i,], paired=TRUE)$statistic)
  }
  
  tab <- data.frame(protein=protnames, pvalue=as.numeric(pvalues), statistic=as.numeric(statistics))
  # set NA values, e.g due to no variance in one of the metrices, to 0
  tab[is.na(tab)] <- 0 
  tab
}

#' @name truncateDecimals 
#' @title Truncate float to defined number of decimal values
#' @description Truncate float to defined number of decimal values
#' @param x float
#' @param digits Number of decimal values
#' @return Truncated number
#' @export
truncateDecimals <- function(x, digits = 2) {
  up <- 10 ^ digits
  x <- x * up
  y <- floor(x)
  y / up
}

#' @name oneSidedTest 
#' @title Recalulate p value from two-sided to one-sided 
#' @description Recalulate p value from two-sided to one-sided 
#' @param sign_value P value from two-sided significance test
#' @param z_value Z value from two-sided significance test
#' @return P value from one sided significance test
#' @export
oneSidedTest <- function(sign_value, z_value){
  if (z_value < 0){
    sign_value <- 1 - (sign_value / 2)
  } else {
    sign_value <- sign_value / 2
  }
  sign_value
}