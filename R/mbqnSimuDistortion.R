#' Perturbation of sample mean and scale
#'
#' @description \code{mbqnSimuDistortion} adds a random perturbation of mean and
#' scale to each column of a matrix.
#' @param x a matrix or data frame.
#' @param s.mean scatter of relative change of mean.
#' @param s.scale scatter of realtive change in scale, i.e. 0.01
#' corresponds to 1 percent.
#' @details Shift and scale the sample mean and standard deviation of a matrix.
#' The perturbation of center and scale relative to mean and standard deviation
#' of each sample are drawn from a Gaussian distribution \eqn{|N(0,\sigma^2)|}
#' with \eqn{\sigma_mean=}\code{s.mean} and
#' \eqn{\sigma_scale}=\code{s.scale}, respectively.
#' @return List with:
#' \item{\code{x.mod}}{perturbed matrix}
#' \item{\code{mx.offset}}{numeric array of shifts of the sample means}
#' \item{\code{mx.scale}}{numeric array of relative scales of the sample
#' standard deviations.}
#' @seealso [mbqnSimuData()] for data generation.
#' @references Brombacher, E., Schad, A., Kreutz, C. (2020). Tail-Robust 
#' Quantile Normalization. BioRxiv.
#' @examples
#' set.seed(1234)
#' x <- mbqnSimuData("omics.dep")
#' df <- mbqnSimuDistortion(x)
#' @author Ariane Schad
#' @export mbqnSimuDistortion
mbqnSimuDistortion <- function(x, s.mean = 0.05, s.scale = 0.01){

    nc <- ncol(x)
    mx.offset <- abs(rnorm(nc))*s.mean+1
    mx.scale <- abs(rnorm(nc))*s.scale+1

    # distort matrix
    mx <- apply(x,2,mean,na.rm =TRUE)
    # x.scaled <- sapply(seq_len(nc), function(i) (x[,i]-mx[i])*mx.scale[i])
    x.scaled <- vapply(seq_len(nc), function(i) (x[,i]-mx[i])*mx.scale[i],
                    FUN.VALUE = numeric(nrow(x)))

    # x.mod <- sapply(seq_len(nc), function(i) (x.scaled[,i]+(mx*mx.offset)[i]))
    x.mod <- vapply(seq_len(nc), function(i) (x.scaled[,i]+(mx*mx.offset)[i]),
                    FUN.VALUE = numeric(nrow(x.scaled)))

    return(distort = list(x.mod = x.mod, mx.offset = mx.offset,
                        mx.scale = mx.scale))
}
