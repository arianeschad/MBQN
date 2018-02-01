#' Mean-balanced quantile normalization
#'
#' Function to compute a mean-balanced quantile normalization
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param FUN mean or median, or an array with dim = numbers of rows of dat. If left empty, quantile normalization is applied without balancing the data
#' @return A matrix of mean- or median-balanced quantile normalized data
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile normalization for processing label-free quantitative proteomics data with abundance-isolated proteins. Biostatistics xxx in prep. 
#' @examples mbqn(dat, mean)
#' mbqn(dat,median)
#' mbqn(dat, x)
#' @description This function uses normalize.quantiles() from the package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
#' @author A. Schad, Aug. 2017


mbqn <- function(dat, FUN = NULL){
  print("Compute mean balanced QN")
  if(!is.null(FUN)){
    print("Compute mean balanced QN")
    library(preprocessCore)
    
    if(is.function(FUN)){
    # sample mean for each row (protein)
    mdat <- apply(dat,1,FUN,na.rm=TRUE)
    }else if(is.function(FUN)){
      mdat <- FUN
      }
    #mean balanced quantile normalisation
    qn_dat <- normalize.quantiles(dat-mdat) + mdat
  }
  else {
  print("QN without mean balancing.")
    library(preprocessCore)
    #quantile normalisation
    qn_dat <- normalize.quantiles(dat)
  }
  return(qn_dat)
}


#mbqn <- function(dat, FUN = mean_fun){
#    mdat <- apply(dat,1,FUN,na.rm=TRUE)
#   qn_dat <- normalize.quantiles(dat-mdat) + mdat
#  return(qn_dat)
#}
