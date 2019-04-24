#' Demonstrate influence of normalization on statistical inference
#'
#' @description Apply a two-sided t-test before and after application
#' of different normalizations to a simulated, differentially expressed and
#' distorted RI feature.
#' @param show.fig Logical that indicates whether results are plot to pdf.
#' @details Apply a two-sided t-test to an RI feature before and after
#' normalization with MBQN. The feature is obtained from a simulated dataset
#' where each sample is distorted in mean and scale.
#' @return Multiple figures of boxplots of intensity distribution for different
#' normalization and of the RI expression profile before and after
#' normalization, t-test results, and a list with:
#' \item{\code{x.mod}}{distorted matrix}
#' \item{\code{mx.offset}}{numeric array of shifts of the sample means}
#' \item{\code{mx.scale}}{numeric array of relative scales of the sample
#' standard deviations}
#' \item{\code{ttest.undistorted}}{output from t-test for the undistorted
#' RI feature}
#' \item{\code{ttest.distorted}}{output from t-test for the distorted
#' RI feature}
#' \item{\code{ttest.mbqndistorted}}{output from t-test for the distorted RI
#' feature after mean balanced quantile normalization.}
#' @seealso [mbqnSimuData()] for data generation and [mbqnSimuDistortion()]
#' for distortion of data.
#' @importFrom utils read.csv untar unzip
#' @importFrom stats rnorm t.test
#' @importFrom graphics image layout points rect matplot
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean/median-balanced
#' quantile normalization. In prep. 2019.
#' @examples
#'\dontrun{
#' set.seed(1234)
#' mbqnDemoTtest(show.fig = TRUE)
#' }
#' @author Ariane Schad
#' @export mbqnDemoTtest
mbqnDemoTtest <- function(show.fig = FALSE){

  mtx <- mbqnSimuData("omics.dep", show.fig = FALSE)
  mtx.mod <- mbqnSimuDistortion(mtx, s.mean = 0.05, s.scale = 0.01)
  bla <- mtx.mod
  mtx.mod <- mtx.mod$x.mod

  res <- mbqnGetNRIfeatures(mtx.mod, low_thr = 0.5, verbose = FALSE)

  # undistorted feature
  feature1 <- mtx[1,]
  # distorted feature
  feature1mod = mtx.mod[1,]
  # feature after normalization
  qn.feature1 = mbqn(mtx.mod, verbose = FALSE)[1,]
  qn.mtx = mbqn(mtx.mod,verbose = FALSE)

  mbqn.mtx = mbqn(mtx.mod, FUN = "mean",verbose = FALSE)
  mbqn.feature1 = mbqn(mtx.mod, FUN = "mean",verbose = FALSE)[1,]

  # Apply t-test:
  # undistorted feature
  ttest.res0 <- t.test(feature1[seq_len(9)], feature1[c(10:18)],
                       var.equal =TRUE)
  # distorted feature
  ttest.res1 <- t.test(feature1mod[seq_len(9)], feature1mod[c(10:18)],
                       var.equal =TRUE)
  # mbqn normalized distorted feature
  ttest.res <- t.test(mbqn.feature1[seq_len(9)], mbqn.feature1[c(10:18)],
                      var.equal =TRUE)

  if(show.fig){
    # compare qn, mbqn and original feature
    par(mfrow = c(1,1))
    mbqnBoxplot(mtx,ylim = c(25,36),
                vals = data.frame(RI = mtx[1,],
                                  NRI = mtx[as.numeric(names(res$nri)[2]),]),
                y.intersp= 1.5)
    dev.copy2pdf(file=file.path(getwd(),"fig_undistorted_boxplot_simu.pdf"),
                 width=8,height=4,paper="a4r",out.type = "pdf")

    #dev.off()
    mbqnBoxplot(mtx.mod,ylim = c(25,39),
                vals = data.frame(RI = mtx.mod[1,],
                                  NRI = mtx.mod[as.numeric(names(res$nri)[2]),]),
                y.intersp= 1.5)
    dev.copy2pdf(file=file.path(getwd(),"fig_distortion_boxplot_simu.pdf"),
                 width=8,height=4,paper="a4r",out.type = "pdf")

    #dev.off()
    mbqnBoxplot(qn.mtx,ylim = c(25,36),
                vals = data.frame(RI = qn.mtx[1,],
                                  NRI = qn.mtx[as.numeric(names(res$nri)[2]),]),
                y.intersp= 1.5)
    dev.copy2pdf(file=file.path(getwd(),"fig_normalized_boxplot_simu.pdf"),
                 width=8,height=4,paper="a4r",out.type = "pdf")

    #dev.off()
    mbqnBoxplot(mbqn.mtx,ylim = c(25,36),
                vals = data.frame(RI = mbqn.mtx[1,],
                                  NRI = mbqn.mtx[as.numeric(names(res$nri)[2]),]),
                y.intersp= 1.5)
    dev.copy2pdf(file=file.path(getwd(),"fig_mbqn_boxplot_simu.pdf"),
                 width=8,height=4,paper="a4r",out.type = "pdf")

    #dev.off()
    matplot(t(rbind(feature1 = feature1,
                    feature1.mod = (feature1mod-mean(feature1mod))/25+mean(feature1),
                    qn.feature1 = (qn.feature1-mean(qn.feature1))+mean(feature1),
                    mbqn.feature1 = (mbqn.feature1-mean(mbqn.feature1))+mean(feature1))),
            type = "b", lty = c(1,1,1), pch = "o", ylab = "intensity", xlab = "sample",
            main = "Differentially expressed RI feature",
            ylim = c(34.48,34.85))
    legend(x=11,y= 34.86, legend = c("feature","distorted feature/25" ,
                                     "QN feature", " MBQN feature"),pch = 1,
           col = c(1,2,3,4), lty= c(1,1,1,1), bty = "n", y.intersp = 1.5,
           x.intersp = 0.2)
    legend(x = .1, y = 34.6,
           legend = paste("p-value (t-test) =",round(ttest.res1$p.value,2),
                          "\np-value (t-test, mbqn) =", round(ttest.res$p.value,4)),
           bty = "n", x.intersp = 0)
    dev.copy2pdf(file=file.path(getwd(),"fig_features_simu.pdf"),
                 width=8,height=4.5,paper="a4r",out.type = "pdf")
  }

  if(ttest.res$p.value<0.05)
    message("H0 (=equal mean) is rejected!")



   return(list(bla$x.mod, bla$mx.offset, bla$mx.scale,
               ttest.undistorted = ttest.res0,
               ttest.distorted = ttest.res1,
               ttest.mbqndistorted = ttest.res))

}
