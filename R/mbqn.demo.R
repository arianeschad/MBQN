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
#' @export
# Installation of package preprocessCore necessary!
# It contains a function for (standard) quantile normalization:
# source('http://bioconductor.org/biocLite.R')
# biocLite('preprocessCore')


mbqn.demo <- function(dat = NULL){

  # if no matrix is given, create a simple dummy matrix
  if(is.null(dat)){
      dat <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3)

    print(dat)
    #     [,1] [,2] [,3]
    #[1,]    5    4    3
    #[2,]    2    1    4
    #[3,]    3    4    6
    #[4,]    4    2    8
  }

  # quantile normalisation
  qn_dat <- preprocessCore::normalize.quantiles(dat)

  # now perform mean balanced qn
  mdat <- mbqn(dat,FUN = median)
  # sample mean for each row (protein)
  mdat <- apply(dat,1,mean,na.rm=TRUE)

  # mean balanced quantile normalisation
  mbqn_dat <- preprocessCore::normalize.quantiles(dat-mdat) + mdat

  return(mbqn_dat)

}

