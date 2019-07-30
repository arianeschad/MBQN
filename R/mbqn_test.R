##########################
#Source: github.com/musto101/wilcox_R

elimna <- function(m){
  #
  # remove any rows of data having missing values
  #
  if(is.null(dim(m)))m<-as.matrix(m)
  ikeep<-c(1:nrow(m))
  for(i in 1:nrow(m))if(sum(is.na(m[i,])>=1))ikeep[i]<-0
  elimna<-m[ikeep[ikeep>=1],]
  elimna
}

olshc4 <- function(x,y,alpha=.05,CN=FALSE,xout=FALSE,outfun=outpro,HC3=FALSE,plotit=FALSE,xlab = "X", ylab = "Y", zlab = "Z",...){
  #
  # Compute confidence intervals via least squares
  # regression using heteroscedastic method
  # recommended by Cribari-Neto (2004).
  # CN=F, degrees of freedom are n-p
  # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
  # All indications are that CN=F is best for general use.
  #
  #  HC3=TRUE, will replace the HC4 estimator with the HC3 estimator.
  #
  x<-as.matrix(x)
  pnum=ncol(x)
  if(nrow(x) != length(y))stop("Length of y does not match number of x values")
  m<-cbind(x,y)
  m<-elimna(m)
  y<-m[,ncol(x)+1]
  x=m[,1:ncol(x)]
  n=length(y)
  nrem=n
  n.keep=length(y)
  x<-as.matrix(x)
  if(xout){
    flag<-outfun(x,plotit=FALSE,...)$keep
    x<-as.matrix(x)
    x<-x[flag,]
    y<-y[flag]
    n.keep=length(y)
    x<-as.matrix(x)
  }
  temp<-lsfit(x,y)
  x<-cbind(rep(1,nrow(x)),x)
  xtx<-solve(t(x)%*%x)
  h<-diag(x%*%xtx%*%t(x))
  n<-length(h)
  d<-(n*h)/sum(h)
  for(i in 1:length(d)){
    d[i]<-min(4, d[i])
  }
  if(HC3)d=2
  hc4<-xtx%*%t(x)%*%diag(temp$res^2/(1-h)^d)%*%x%*%xtx
  df<-nrow(x)-ncol(x)
  crit<-qt(1-alpha/2,df)
  if(CN)crit=qnorm(1-alpha/2)
  al<-ncol(x)
  p=al-1
  ci<-matrix(NA,nrow=al,ncol=7)
  lab.out=rep("Slope",p)
  dimnames(ci)<-list(c("(Intercept)",lab.out),c("Coef.","Estimates",
                                                "ci.lower","ci.upper", "t-value", "p-value","Std.Error"))
  for(j in 1:al){
    ci[j,1]<-j-1
    ci[j,2]<-temp$coef[j]
    ci[j,3]<-temp$coef[j]-crit*sqrt(hc4[j,j])
    ci[j,4]<-temp$coef[j]+crit*sqrt(hc4[j,j])
    test<-temp$coef[j]/sqrt(hc4[j,j])
    ci[j,5]<- test
    ci[j,6]<-2*(1-pt(abs(test),df)) # pt: Student t distribution, lower.tail = TRUE
    if(CN)ci[j,6]<-2*(1-pnorm(abs(test),df))
  }
  ci[,7]=sqrt(diag(hc4))
  if(plotit){
    if(pnum==1){
      plot(x[,-1],y,xlab=xlab,ylab=ylab)
      abline(ci[,2])
    }
    if(pnum==2){
      regp2plot(x[,-1],y,regfun=ols,xlab=xlab,ylab=ylab,zlab=zlab)
    }}
  list(n=nrem,n.keep=n.keep,ci=ci, cov=hc4)
}


pcorhc4 <- function(x,y,alpha=.05,CN=FALSE){
  #
  #   Compute a .95 confidence interval for Pearson's correlation coefficient.
  #   using the HC4 method
  #
  # CN=F, degrees of freedom are n-p; seems better for general use.
  # CN=T  degrees of freedom are infinite, as done by Cribari-Neto (2004)
  #
  xy<-elimna(cbind(x,y))
  x<-xy[,1]
  y<-xy[,2]
  z1=(x-mean(x))/sqrt(var(x))
  z2=(y-mean(y))/sqrt(var(y))
  ans=olshc4(z1,z2,alpha=alpha,CN=CN)
  list(r=cor(x,y),ci=ans$ci[2,3:4],t.statistic=ans$ci[2,5], p.value=ans$ci[2,6])
}

comdvar <- function(x,y,alpha=.05){
  #
  # Test the hypothesis that two dependent variables have equal variances.
  # A heteroscedastic version of the Morgan-Pitman test is used.
  # (The HC4 estimator is used to deal with heteroscedasticity)
  #
  # Returns: p-value, variance of x, variance of y
  #
  xy=elimna(cbind(x,y))
  est1=var(xy[,1])
  est2=var(xy[,2])
  pv=pcorhc4(xy[,1]-xy[,2],xy[,1]+xy[,2],alpha=alpha)
  list(p.value=pv$p.value, statistic=pv$t.statistic, est1=est1, est2=est2)
}

#################################
library("SummarizedExperiment")
getMtx <- function(pxd_id, meanMedian){
  # Load file
  out <- mbqnLoadFile(pxd_id, file.pattern = "proteinGroups.txt");
  
  # filter for potential contaminants and identified only by site features
  out <- out[!rowData(out)[["ixs"]],] 
  
  # extract data and feature annotation
  mtx <- assays(out)[["data"]]
  
  colnames(mtx) <- gsub("LFQ intensity","",colnames(mtx))
  #mtx[is.na(mtx)] <- 0 # added EB
  
  # remove all proteins with NAs
  mtx <- na.omit(mtx)
  
  mbqn.mtx <- mbqn(mtx,FUN = meanMedian)
  qn.mtx <- mbqn(mtx,FUN = NULL)
  
  row.names(mbqn.mtx) <- row.names(mtx)
  row.names(qn.mtx) <- row.names(mtx)
  list(mtx=mtx, mbqn.mtx=mbqn.mtx, qn.mtx=qn.mtx)
}


library(PairedData)
getPvalue <- function(mtx1, mtx2, method){
  protnames <- c()
  pvalues <- c()
  statistics <- c()
  options("scipen"=100, "digits"=4)
  
  if (is.null(row.names(mtx1))){
    rownames(mtx1) <- 1:nrow(mtx1)
  }
  
  for (i in seq_len(nrow(mtx))){
    protnames <- c(protnames, row.names(mtx1)[i])
    if (method == "pitman"){
      pvalues <- c(pvalues, Var.test(mtx1[i,], mtx2[i,], paired=TRUE)$p.value)
      statistics <- c(statistics, Var.test(mtx1[i,], mtx2[i,], paired=TRUE)$statistic)
    } else if (method == "hc4"){
      pvalues <- c(pvalues, comdvar(mtx1[i,], mtx2[i,])$p.value)
      statistics <- c(statistics, Var.test(mtx1[i,], mtx2[i,])$statistic)
    }
  }
  tab <- data.frame(protein=protnames, pvalue=as.numeric(pvalues), statistic=as.numeric(statistics))
  # set NA values, e.g due to no valiance in one of the metarices, to ß
  tab[is.na(tab)] <- 0 
  tab
}

logit_norm <- function(df_qn, df_mbqn, threshold, pxd_id2, meanMedian2){
  set.seed(2)
  
  df_qn$sign <- ifelse(df_qn$pvalue < thr, 1, 0)
  df_mbqn$sign <- ifelse(df_mbqn$pvalue < thr, 1, 0)
  
  model_qn <- glm(sign ~ nri_freq, data = df_qn, family = "binomial")
  model_mbqn <- glm(sign ~ nri_freq, data = df_mbqn, family = "binomial")
  
  print(summary(model_qn))
  z_value_qn <- summary(model_qn)$coefficients[2,3]
  sign_value_qn <- summary(model_qn)$coefficients[2,4]
  
  if (z_value_qn < 0){
    sign_value_qn <- 1 - (sign_value_qn / 2)
  } else {
    sign_value_qn <- sign_value_qn / 2
  } 
  
  print(paste("Sign right-tailed QN: ", round(sign_value_qn, 3)))
  
  print(summary(model_mbqn))
  z_value_mbqn <- summary(model_mbqn)$coefficients[2,3]
  sign_value_mbqn <- summary(model_mbqn)$coefficients[2,4]
  
  if (z_value_mbqn < 0){
    sign_value_mbqn <- 1 - (sign_value_mbqn / 2)
  } else {
    sign_value_mbqn <- sign_value_mbqn / 2
  }
  
  print(paste("Sign right-tailed MBQN: ", round(sign_value_mbqn, 3)))
  
  predProbs_qn = predict(model_qn, type="response") ## extract fitted probabilities
  predProbs_mbqn = predict(model_mbqn, type="response") ## extract fitted probabilities
  
  f1 <- approxfun(predProbs_qn - predProbs_mbqn, df_qn$nri_freq, rule=2)
  intersect <- f1(0)

  if (sign_value_qn < 0.01) { #0.001
    if (intersect == min(df_qn$nri_freq)){
      warning('QN Problem: Use MBQN for whole dataset.')
    } else {
      warning(paste('QN Problem: Use MBQN. NRI-threshold =',intersect,'suggested.'))
    }
  }
  
  theme_set(theme_minimal())
  plot_combined <- ggplot() +
    geom_point(aes(x=df_qn$nri_freq, df_qn$sign), color="#CC6666", alpha=0.3) +
    geom_line(aes(x=df_qn$nri_freq, y=predProbs_qn), color="#CC6666") +
    geom_point(aes(x=df_mbqn$nri_freq, df_mbqn$sign), color="#66CC99", alpha=0.3) +
    geom_line(aes(x=df_mbqn$nri_freq, y=predProbs_mbqn), color="#66CC99") +
    ylim(0,1) +
    labs(
      x = "Frequency of RI proteins [%]",
      y = "Probability of significant p-values",
      title = paste(pxd_id2, meanMedian2, "Threshold:", threshold, ", QN: red (", round(sign_value_qn, 3),"), MBQN: green (", round(sign_value_mbqn, 3), ")"))
  ggsave(paste0("logistRegr_nriVsProbab_QNandMBQN_thresh_", threshold,"_", pxd_id2, "_", meanMedian2, ".pdf"), plot = plot_combined, width=10, height=8)
}




meanMedian <- "median" # mean, median
pxd_id <- "PXD001584"  
#pxd_id <- "PXD005861"
#pxd_id <- "PXD006617"
##### pxd_id <- "PXD005138"
#pxd_id <- "simulated"
mtx <- mbqnSimu(offset=100)
mbqn.mtx <- mbqnNRI(x = mtx, FUN = median, verbose = FALSE) # mbqn
qn.mtx <- mbqnNRI(x = mtx, FUN = NULL, verbose = FALSE) # QN


library(ggplot2)
library(reshape2)

# list <- getMtx(pxd_id, meanMedian)
# 
# mtx <- list$mtx
# mbqn.mtx <- list$mbqn.mtx
# qn.mtx <- list$qn.mtx

pvalues_mbqn_hc4 <- getPvalue(mtx, mbqn.mtx, "hc4")
pvalues_qn_hc4 <- getPvalue(mtx, qn.mtx, "hc4")

low_thr <- 0.0
res  <- mbqnGetNRIfeatures(mtx, low_thr = low_thr)
nri_df <- data.frame(res$nri)
nri_freq <- c()
for (row in seq_len(nrow(nri_df))){
  nri <- nri_df[row,]$Var1
  freq <- nri_df[row,]$Freq
  nri_freq <- c(nri_freq, freq)
  prot <- row.names(mtx)[nri]
}

combined_mbqn_hc4 <- data.frame(nri_freq, pvalues_mbqn_hc4)
combined_qn_hc4 <- data.frame(nri_freq, pvalues_qn_hc4)

thr = 0.01
logit_norm(combined_qn_hc4, combined_mbqn_hc4, thr, pxd_id, meanMedian)
