#' Get the k largest/smallest elements
#'
#' @description Extract the k largest or smallest values and their indices for
#' each column of a matrix.
#' @param x a data matrix or data frame.
#' @param k an integer specifying the number of extreme values. Must be
#' \code{<= nrows(x)}.
#' @param flag use "min" or "max" (default) to select smallest or largest
#' elements.
#' @details Order the values of each column of \code{x} and determine the
#' k smallest (\code{flag = "min"}) or largest (\code{flag = "max"}) values and
#' their indices. NA's in the data are ignored.
#' @return List with elements:
#' \item{\code{ik}}{indices of ordered extreme values}
#' \item{\code{minmax}}{ordered extreme values.}
# #' @concept quantile, quantile normalization, rank invariance
#' @references Schad, A. and Kreutz, C., MBQN: R package for
#' mean/median-balanced quantile
#' normalization. In prep. 2019
#' @examples
#' x <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3) # Create a data matrix
#' getKminmax(x, k = 5, "max") # get indices of the 5 largest values in each column
#' @author Ariane Schad
#  Aug. 2017
#' @export getKminmax
getKminmax <- function(x,k,flag = "max"){

  if(flag == "min"){
    decreasing <- FALSE
  }else if(flag=="max"){
    decreasing <- TRUE
  }

  ## sort ##
  na.last <- TRUE
  ik <- apply(x,2,order, decreasing = decreasing, na.last = na.last )
  minmax <- apply(x,2,function(y)y[order(y, decreasing = decreasing,
                                         na.last = na.last )])
  ## select the first k elements ##
  ik <- ik[seq_len(k),]
  minmax <- minmax[seq_len(k),]

  ## return results
  return(list(ik = ik, minmax = minmax))
}
