#' Simulate data matrix as input for the mean-balanced quantile normalization
#'
#' @description Generate a data matrix for illustration of mean-balanced quantile normalization
#' @param N number of rows of data matrix
#' @param M number of columns of data matrix
#' @return N x M matrix 
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile normalization for processing label-free quantitative proteomics data with abundance-isolated proteins. Biostatistics xxx in prep. 
#' @examples 
#' mbqn.simu(1000,10)
#' @details Matrix entries are from a standard normal distribution of mean 0 and sigma 1. The data is used for mbqn.demo()
#' @author A. Schad, Aug. 2017

# begin function
mbqn.simu_dat <- function(N = NULL, M = NULL){
  
  if(is.null(N)) N <- 1000
  if(is.null(M)) M <- 10
  
  sample.size <- M
  mu <- rep(0,M)
  sigma <- rep(1,M)
  set.seed(123)
  dat <- mapply(function(x,y){rnorm(x,y,n=N)},x=mu,y=sigma)
  dim(dat)
  return(dat)
  }




