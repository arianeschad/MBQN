#' Demonstration of mean/median-balanced quantile normalization
#'
#' @description This function demonstrates mean/median-balanced quantile normalization
#' @param dat a matrix where rows represent features, e.g. protein
#' abundances/intensities, and columns represent experimental samples.
#' Default \code{NULL} - a simple matrix is generated.
# where dat is a data matrix where columns correspond
# to experimental samples or (technical/biological) replicates and
# rows correspond to features/proteins/peptides
#' @return Various graphics illustrating saturation in case of
#' quantile normalization of high expression for a single protein
#' abundances across samples
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile
#' normalization for processing label-free quantitative proteomics
#' data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples
#' mbqnDemo()
#' dat <- mbqnSimuData()
#' mbqnDemo(dat)
#' @details Normalize a matrix and return boxplots of quantile normalized and mean balanced normalized
#' data. An omics-like data matrix can be generated with \code{mbqn.simu.dat()}.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
# Aug. 2017
#' @export mbqnDemo

mbqnDemo <- function(dat = NULL){

  # if no matrix is given, create a simple dummy matrix
  if(is.null(dat)){
    dat <- matrix(c(5,2,3,NA,2,4,1,4,2,3,1,4,6,NA,1,3,NA,1,4,3,NA,1,2,3),ncol=4)

    print(dat)
    #      [,1] [,2] [,3] [,4]
    # [1,]    5    1    6    4
    # [2,]    2    4   NA    3
    # [3,]    3    2    1   NA
    # [4,]   NA    3    3    1
    # [5,]    2    1   NA    2
    # [5,]    4    4    1    3

  }

  # perform qn, median balanced qn, and qn with median balanced nri feature
  qn_dat <- mbqn(dat,FUN=NULL)
  mbqn_dat <- mbqn(dat,FUN = median)
  qn_nri_dat <- mbqnNRI(dat,FUN = "median", low_thr = 0.5)

  # sample mean for each row (protein)
  #mdat <- apply(dat,1,mean,na.rm=TRUE)

  # check saturation i.e. for rank invariance
  res <- mbqnCheckSaturation(dat, save_fig = FALSE, verbose = FALSE)
  ylim <- range(dat, na.rm =T)

  plot.new()
  frame()
  par(mfrow=c(2,2))
  # create a boxplot for dat
  mbqnBoxplot(dat, filename = NULL, add.leg = F, main = "data", ylim = ylim)
  # create a boxplot for qn-data
  mbqnBoxplot(qn_dat, filename = NULL, irow = res$ip, add.leg = F, main = "QN data", ylim = ylim)
  # create a boxplot for mbqn-data
  mbqnBoxplot(mbqn_dat, filename = NULL, add.leg = F, main = "MBQN data", ylim = ylim)
  # create a boxplot for qn-data with nri features median balanced
  mbqnBoxplot(qn_nri_dat, irow = res$ip, filename = NULL, add.leg = F, main ="QN with MBQN of NRI", ylim = ylim)

}

