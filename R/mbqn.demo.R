#' Demonstration of mean/median-balanced quantile normalization
#'
#' @description This function demonstrates mean-balanced quantile normalization
#' @param dat A data matrix where rows represent features, e.g. protein
#' abundances/intensities, and columns experimental samples or
#' (technical/biological) replicates.
#' Default NULL - a simple matrix is generated.
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
#' mbqn.demo
#' mbqn.demo(dat)
#' @details This function uses \code{normalize.quantiles()} from the package
#' preprocessCore that can be installed from http://bioconductor.org/biocLite.R by
#' source('http://bioconductor.org/biocLite.R')
#' biocLite('preprocessCore'). A data matrix can be computed with mbqn.simu_dat().
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' Aug. 2017
#' @export mbqn.demo
# Installation of package preprocessCore necessary!
# It contains a function for (standard) quantile normalization:
# source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore')


mbqn.demo <- function(dat = NULL){

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
  qn_nri_dat <- mbqn.nri(dat,FUN = median, low_thr = 0.6)

  # sample mean for each row (protein)
  #mdat <- apply(dat,1,mean,na.rm=TRUE)

  # check saturation i.e. for rank invariance
  res <- mbqn.check_saturation(dat)

  plot.new()
  frame()
  # create a boxplot for dat
  mbqn.boxplot(dat)
  # create a boxplot for qn-data
  mbqn.boxplot(qn_dat)
  # create a boxplot for mbqn-data
  mbqn.boxplot(mbqn_dat)
  # create a boxplot for qn-data with nri features median balanced
  mbqn.boxplot(qn_nri_dat, irow = res$ip)

  return(mbqn_dat)

}

