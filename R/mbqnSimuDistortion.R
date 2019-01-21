#' Simulate distortion of mean and scale of a data matrix
#'
#' @description Simulate an omics-like, columnwise distortion in mean and scale
#' of a data matrix
#' @param ncol number of columns of data matrix
#' @param seed Seed for random number generator. Default \code{seed <- 1234}
#' @param max_mean Absolute maximum shift of mean from 0 to max_mean
#' @param max_scale Realtive change in scale, i.e. 0.01 corresponds to 1\%
#' @details For model "rand" the random matrix entries are drawn from a standard normal distribution
#' of mean 0 and sigma 1. For model "omics" the meand and standard deviation of each row i is
#' drawn from Gaussian distribution \eqn{N(\mu_i,sd_i^2)} with \eqn{sd_i~N(0,0.0625)} and \eqn{\mu_i~N(28,4)}.
#' About 21\% of the values are omitted according to the missing value pattern present in the LFQ intensities
#' of PXD001584 (Ramond et al. 2015).
#' @return \code{matrix} of size nrow x ncol
#' @importFrom utils read.csv untar unzip
#' @importFrom stats rnorm t.test
#' @importFrom graphics image layout points rect matplot
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019. \cr
#' Ramond, E. et al. (2015) Importance of host cell arginine uptake in Francisella phagosomal
#' escape and ribosomal protein amounts. Mol Cell Proteomics 14, 870-881.
#' @examples
#'\dontrun{
#' MBQN:::mbqnSimuDistortion(ncol = 9)
#' }
#' @author A. Schad, \email{ariane.schad@zbsa.de}
# #' @export mbqnSimuDistortion
mbqnSimuDistortion <- function(ncol = NULL,max_mean = 0.05, max_scale = 0.01, seed = 1234){

  if(!is.null(seed)) set.seed(seed)

  mx.offset <- abs(rnorm(ncol))*max_mean+1
  mx.scale <- abs(rnorm(ncol))*max_scale+1

  # sapply(mtx, function(i) mtx[,i]*mx.offset[i])
  mtx <- mbqnSimuData("omics")
  # add extra feature, that shows diff. expr.
  #mtx <- mbqn(mtx)
  set.seed(1111)
  mtx <- rbind(c(rep(34.7,9),rep(34.6, 9))+rnorm(ncol)*0.02,mtx)

  mx <- apply(mtx,2,mean,na.rm =T)
  mtx.scaled <- sapply(1:ncol, function(i) (mtx[,i]-mx[i])*mx.scale[i])
  mtx.mod <- sapply(1:ncol, function(i) (mtx.scaled[,i]+(mx*mx.offset)[i]))

  res <- mbqnCheckSaturation(mtx.mod,FUN= "mean", low_thr = 0.5, verbose = F, show_nri_only = T, y.intersep = 0.8, save_fig = T)

  # selected feature
  feature1 = mtx[1,]
  feature1mod = mtx.mod[1,]
  qn.feature1 = mbqn(mtx.mod, verbose = F)[1,]
  qn.mtx = mbqn(mtx.mod,verbose = F)

  mbqn.mtx = mbqn(mtx.mod, FUN = "mean",verbose = F)
  mbqn.feature1 = mbqn(mtx.mod, FUN = "mean",verbose = F)[1,]

  # apply t-test
  ttest.res <- t.test(mbqn.feature1[1:9], mbqn.feature1[10:18],var.equal =TRUE)
  ttest.res0 <- t.test(feature1[1:9], feature1[10:18],var.equal =TRUE)$p.value
  ttest.res1 <- t.test(feature1mod[1:9], feature1mod[10:18],var.equal =TRUE)$p.value

  # compare qn, mbqn and original feature
  dev.off()
  mbqnBoxplot(mtx,ylim = c(22,36), vals = data.frame(RI = mtx[1,],NRI = mtx[as.numeric(names(res$nri)[2]),]), y.intersp= 0.8)
  dev.copy2pdf(file=file.path(getwd(),"fig_undistorted_boxplot_simu.pdf"),width=8,height=4,paper="a4r",out.type = "pdf")

  dev.off()
  mbqnBoxplot(mtx.mod,ylim = c(21,39), vals = data.frame(RI = mtx.mod[1,],NRI = mtx.mod[as.numeric(names(res$nri)[2]),]), y.intersp= 0.8)
  dev.copy2pdf(file=file.path(getwd(),"fig_distortion_boxplot_simu.pdf"),width=8,height=4,paper="a4r",out.type = "pdf")

  dev.off()
  mbqnBoxplot(qn.mtx,ylim = c(22,36), vals = data.frame(RI = qn.mtx[1,],NRI = qn.mtx[as.numeric(names(res$nri)[2]),]), y.intersp= 0.8)
  dev.copy2pdf(file=file.path(getwd(),"fig_normalized_boxplot_simu.pdf"),width=8,height=4,paper="a4r",out.type = "pdf")

  dev.off()
  mbqnBoxplot(mbqn.mtx,ylim = c(22,36), vals = data.frame(RI = mbqn.mtx[1,],NRI = mbqn.mtx[as.numeric(names(res$nri)[2]),]), y.intersp= 0.8)
  dev.copy2pdf(file=file.path(getwd(),"fig_mbqn_boxplot_simu.pdf"),width=8,height=4,paper="a4r",out.type = "pdf")

  dev.off() #-mean(feature1mod))/50+mean(qn.feature1)
  matplot(t(rbind(feature1 = feature1,
                  feature1.mod = (feature1mod-mean(feature1mod))/50+mean(feature1),
                  qn.feature1 = (qn.feature1-mean(qn.feature1))+mean(feature1),
                  mbqn.feature1 = (mbqn.feature1-mean(mbqn.feature1))+mean(feature1))),
          type = "b", lty = c(1,1,1), pch = "o", ylab = "intensity", xlab = "sample",
          main = "Differentially expressed feature",
          ylim = c(34.48,34.85))
  legend(x=11,y= 34.88, legend = c("feature","distorted feature/50" ,"QN feature", " MBQN feature"),pch = 1,
         col = c(1,2,3,4), lty= c(1,1,1,1), bty = "n", y.intersp = 0.8, x.intersp = 0.2)
  legend(x = -.5, y = 34.6,
         legend = paste("p-value (t-test) =",round(ttest.res1,2), "\np-value (t-test, mbqn) =", round(ttest.res$p.value,4)),
         bty = "n", x.intersp = 0)
  dev.copy2pdf(file=file.path(getwd(),"fig_features_simu.pdf"),width=8,height=4.5,paper="a4r",out.type = "pdf")


  if(ttest.res$p.value<0.05)
    print("H0 (=equal mean) is rejected!")



  return(distort = list(mx.offset = mx.offset, mx.scale = mx.scale))

}
