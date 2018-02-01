#' Demo of mean/median-balanced quantile normalization
#'
#' @description This function demonstrates mean-balanced quantile normalization
#' @param dat Default NULL. Optionally a user-supplied data matrix with rows - features, e.g. protein abundances; columns - samples. Data matrices can be computed with
#' mbqn.simu_dat()
#' @return Various graphics illustrating saturation in case of quantile normalization of high expression for a single protein abundances across samples
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile normalization for processing label-free quantitative proteomics data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples
#' demo_mbqn
#' demo_mbqn(dat) where dat is a data matrix where columns correspond to experimental samples or (technical/biological) replicates and rows correspond to features/proteins/peptides
#' @details This function uses normalize.quantiles() from the package preprocessCore that can be installed from http://bioconductor.org/biocLite.R by
#' source('http://bioconductor.org/biocLite.R')
#' biocLite('preprocessCore')
#' @author A. Schad, Aug. 2017


#install the package, it contains a function for quantile normalization:
#source('http://bioconductor.org/biocLite.R')
#biocLite('preprocessCore')

mbqn.demo <- function(dat = NULL){

  library(preprocessCore)
  if(is.null(dat)){
    # the function expects a matrix as input
    # create a matrix using the same example
    # each column corresponds to a sample
    # each row corresponds to a protein

    dat <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),
                  ncol=3)
    dat
    #     [,1] [,2] [,3]
    #[1,]    5    4    3
    #[2,]    2    1    4
    #[3,]    3    4    6
    #[4,]    4    2    8
  }
  dim(dat)

  #quantile normalisation
  qn_dat <- normalize.quantiles(dat)

  # now perform mean balanced qn

  # sample mean for each row (protein)
  mdat <- apply(dat,1,mean,na.rm=TRUE)

  #mean balanced quantile normalisation
  mbqn_dat <- normalize.quantiles(dat-mdat) + mdat
}

