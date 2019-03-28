#' Generate a random/structured data matrix
#'
#' @description Generate a random data matrix with or without proteomics, log-transformed feature
#' intensity-like properties.
#' @param model character indicating one of the three different type of models: \code{"rand"}(default) a Gaussian random matrix
#' of size nrow x ncol (default 1000 x 10), \code{"omics"} a Gaussian random matrix of size 1264 x 18 that mimics intensity profiles and
#' missing values as present in real data, and \code{"omics.dep"} is the same as \code{"omics"} but with
#' an additional single, differentially expressed RI feature.
#' @param nrow number of rows of data matrix (only for \code{model = "rand"}).
#' @param ncol number of columns of data matrix (only for \code{model = "rand"}).
#' @param seed seed for random number generator.
#' @param show.fig logical inidicating whether data properties are plot to figure (only for \code{model = "omics"} and \code{model = "omics.dep"}).
#' @details For model \code{"rand"}, each matrix element is drawn from a standard normal
#' distribution \eqn{N(0,1)}. For model \code{"omics"}, the matrix elements of each row
#' are drawn from a Gaussian distribution \eqn{N(\mu_i,\sigma_i^2)} where the
#' mean and standard deviation itself are drawn Gaussian distributions, i.e.
#' \eqn{\sigma_i~N(0,0.0625)} and \eqn{\mu_i~N(28,4)}. About 35\% of the matrix
#' values are set to NA according to the missing value pattern present in real protein LFQ
#' intensities. For model \code{"omics.dep"}, a single differentially epxressed
#' RI feature is stacked on top of the matrix from model \code{"omics"}.
#' @return \code{matrix} of size nrow x ncol.
#' @seealso [example_NApattern()] for description of missing value pattern.
#' @importFrom utils read.csv untar unzip
#' @importFrom stats rnorm
#' @importFrom graphics image layout points rect
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean/median-balanced quantile normalization.
#' In prep. 2019
#' @examples
#' mbqnSimuData(model = "rand")
#' mbqnSimuData(model = "rand", 2000,6)
#' mbqnSimuData(model = "omics")
#' mbqnSimuData(model = "omics.dep")
#' @author Ariane Schad
#' @export mbqnSimuData
mbqnSimuData <- function(model = "rand", nrow = NULL, ncol = NULL,seed = 1234, show.fig = FALSE){

  if(is.null(nrow)) nrow <- 1000
  if(is.null(ncol)) ncol <- 10
  if(!is.null(seed)) set.seed(seed)

  if(model=="rand"){
    # generate a random matrix without NAs
    dat <- replicate(ncol, rnorm(nrow))
    if(show.fig){
      image(t(dat), xlab = "sample", ylab = "feature row", main = "simulated data", axes= FALSE)
      axis(1, at = seq(1, ncol(dat), by = 1)/ncol(dat), labels = c(1:ncol(dat)))
    }
  }

  if(model=="omics" || model == "omics.dep"){
    # generate a structured random matrix with NAs
    # the NA structure is extracted from a real dataset

    # load an internal stored MV pattern extracted from a real proteomics dataset from PRIDE
    mtx <- example_NApattern
    mtx[mtx==0] <- NA

    dat <- replicate(dim(mtx)[2], rnorm(dim(mtx)[1])*0.25)+rnorm(dim(mtx)[1])*2+28
    s.dat <-sort(apply(dat, 1,mean, na.rm =TRUE), index.return = TRUE , decreasing = TRUE)
    dat <- dat[s.dat$ix,]
    dat[is.na(mtx)] <- NA

    if(show.fig){
      par(mfrow=c(2,2))
      image(t(mtx), xlab = "sample", ylab = "feature row (sorted)", main = "MV pattern", axes = FALSE)
      axis(1, at = seq(1, ncol(dat), by = 1)/ncol(dat), labels = c(1:ncol(dat)))
      image(t(dat), xlab = "sample", main = "simulated", axes= FALSE)
      axis(1, at = seq(1, ncol(dat), by = 1)/ncol(dat), labels = c(1:ncol(dat)))

      plot(apply(dat,1,mean,na.rm =TRUE),apply(is.na(dat),1,mean,na.rm =TRUE),
           col =1,
           xlab = "mean feature intensity", ylab = "MV frequency")
      legend("topright",legend = c("simulated"),  pch = 1, col = c(1), bty ="n", y.intersp = 1.5)
    }

    if(model=="omics.dep"){
      # add extra feature with diff. expr.
      set.seed(1111)
      ncol <- dim(dat)[2]
      dat <- rbind(c(rep(34.7,9),rep(34.6, 9))+rnorm(ncol)*0.02,dat)

    }

  }

  return(dat)

}




