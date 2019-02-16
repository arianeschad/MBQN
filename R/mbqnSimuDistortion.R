#' Perturbation of sample mean and scale
#'
#' @description \code{mbqnSimuDistortion} adds a random perturbation of mean and scale to each column of a matrix.
#' @param x a matrix or data frame.
#' @param s.mean scatter of relative change of mean.
#' @param s.scale scatter of realtive change in scale, i.e. 0.01 corresponds to 1\%.
#' @param seed seed for random number generator.
#' @details Shift and scale the sample mean and standard deviation of a matrix. The perturbation of center and scale relative to mean and standard deviation of each sample
#' are drawn from a Gaussian distribution \eqn{|N(0,\sigma^2)|} with \eqn{\sigma_mean=}\code{s.mean} and
#' \eqn{\sigma_scale}=\code{s.scale}, respectively.
#' @return List with:
#' \item{\code{x.mod}}{perturbed matrix}
#' \item{\code{mx.offset}}{numeric array of shifts of the sample means}
#' \item{\code{mx.scale}}{numeric array of relative scales of the sample standard deviations.}
#' @seealso [mbqnSimuData()] for data generation.
#' @family data
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples
#'\dontrun{
#' x <- mbqnSimuData("omics.dep")
#' df <- mbqnSimuDistortion(x)
#' }
#' @author Ariane Schad
#' @export mbqnSimuDistortion
mbqnSimuDistortion <- function(x, s.mean = 0.05, s.scale = 0.01, seed = 1234){

  if(!is.null(seed)) set.seed(seed)

  ncol <- dim(x)[2]
  mx.offset <- abs(rnorm(ncol))*s.mean+1
  mx.scale <- abs(rnorm(ncol))*s.scale+1

  # distort matrix
  mx <- apply(x,2,mean,na.rm =T)
  x.scaled <- sapply(1:ncol, function(i) (x[,i]-mx[i])*mx.scale[i])
  x.mod <- sapply(1:ncol, function(i) (x.scaled[,i]+(mx*mx.offset)[i]))

  return(distort = list(x.mod = x.mod, mx.offset = mx.offset, mx.scale = mx.scale))

}
