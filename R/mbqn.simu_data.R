#' Simulate data matrix for mean-balanced quantile normalization
#'
#' @description Generate a data matrix for illustration of mean-balanced.
#' quantile normalization.
#' @param N number of rows of data matrix.
#' @param M number of columns of data matrix.
#' @param seed Seed for random number generator. Default \code{seed <- 1234}.
#' @details Matrix entries are from a standard normal distribution
#' of mean 0 and sigma 1. This function is used for mbqn.demo().
#' @return \code{matrix} of size N x M
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced
#' quantile normalization for processing label-free quantitative
#' proteomics data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples
#' mbqn.simu(1000,10)
#' @author A. Schad, \email{ariane.schad@zbsa.de}

# begin function
mbqn.simu_dat <- function(N = NULL, M = NULL,seed = NULL){

  if(is.null(N)) N <- 1000
  if(is.null(M)) M <- 10
  if(is.null(seed)) seed <- 1234
  sample.size <- M
  mu <- rep(0,M)
  sigma <- rep(1,M)
  set.seed(seed)
  dat <- mapply(function(x,y){rnorm(x,y,n=N)},x=mu,y=sigma)
  dim(dat)
  return(dat)
  }




