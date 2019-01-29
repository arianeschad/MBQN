#' Generate a random/structured data matrix
#'
#' @description Generate a random data matrix with or without proteomics, log-transformed feature
#' intensity-like properties
#' @param model Two different type of matrix models are available: "rand" generates a random matrix
#' of size nrow x ncol, "omics" generates a Gaussian random matrix which mimics intensity profiles and
#' missing values as present in a real data set.
#' @param nrow number of rows of data matrix.
#' @param ncol number of columns of data matrix.
#' @param seed Seed for random number generator. Default \code{seed <- 1234}.
#' @details For model "rand" the matrix entries are drawn from a standard normal
#' distribution \eqn{N(0,1)}. For model "omics" the matrix entries of each row
#' are drawn from a Gaussian distribution \eqn{N(\mu_i,\sigma_i^2)} where the
#' mean and standard deviation itself are drawn Gaussian distributions, i.e.
#' \eqn{\sigma_i~N(0,0.0625)} and \eqn{\mu_i~N(28,4)}. About 35\% of the matrix
#' values are set to NA according to the missing value pattern present in the protein LFQ
#' intensities of PXD001584 \[1\].
#' @return \code{matrix} of size nrow x ncol
#' @family data
#' @seealso [MBQN::example_NApattern()] for description of missing value pattern
#' @importFrom utils read.csv untar unzip
#' @importFrom stats rnorm
#' @importFrom graphics image layout points rect
# #' @keywords quantile normalization proteomics
#' @references
#' \[1\] Ramond, E. et al. (2015) Importance of host cell arginine uptake in Francisella phagosomal
#' escape and ribosomal protein amounts. Mol Cell Proteomics 14, 870-881.\cr
#' @examples
#' \dontrun{
#' mbqnSimuData(model = "rand", 1000,10)
#' mbqnSimuData(model = "rand")
#' mbqnSimuData(model = "omics")
#' }
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export mbqnSimuData
mbqnSimuData <- function(model = "rand", nrow = NULL, ncol = NULL,seed = 1234){

  if(is.null(nrow)) nrow <- 1000
  if(is.null(ncol)) ncol <- 10
  if(!is.null(seed)) set.seed(seed)

  if(model=="rand"){
    # generate a random matrix without NAs
    dat <- replicate(ncol, rnorm(nrow))
  }

  if(model=="omics"){
    # generate a structured random matrix with NAs
    # the NA structure is extracted from a real dataset

    # load an internal stored MV pattern extracted from a real proteomics dataset from PRIDE
    mtx <- MBQN:::example_NApattern
    mtx[mtx==0] <- NA

    dat <- replicate(dim(mtx)[2], rnorm(dim(mtx)[1])*0.25)+rnorm(dim(mtx)[1])*2+28
    s.dat <-sort(apply(dat, 1,mean, na.rm =T), index.return = TRUE , decreasing = TRUE)
    dat <- dat[s.dat$ix,]
    dat[is.na(mtx)] <- NA

    plot.new()
    frame()
    par(mfrow=c(2,2))
    image(t(mtx), xlab = "sample", ylab = "feature row (sorted)", main = "MV pattern", axes = FALSE)
    axis(1, at = seq(1, ncol(dat), by = 1)/ncol(dat), labels = c(1:ncol(dat)))
    image(t(dat), xlab = "sample", main = "simulated", axes= FALSE)
    axis(1, at = seq(1, ncol(dat), by = 1)/ncol(dat), labels = c(1:ncol(dat)))

    plot(apply(dat,1,mean,na.rm =T),apply(is.na(dat),1,mean,na.rm =T),
         col =1,
         xlab = "mean feature intensity", ylab = "MV frequency")
    legend("topright",legend = c("simulated"),  pch = 1, col = c(1), bty ="n", y.intersp = 1.5)
  }

  return(dat)

}




