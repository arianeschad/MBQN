#' Get indices and values of the k largest/smallest elements
#'
#' @description Extract the k largest or smallest values and their indices of each column of a matrix.
#' @param X a data matrix.
#' @param k integer specifying the number of selected extreme values. Must be less \code{<= nrows(X)}.
#' @param flag "min" or "max"(default)
#' @details Order the values of each column of X and determine the
#' k smallest (\code{flag = "min"}) or largest (\code{flag = "max"}) values and their indices. NA's in the
#' data are ignored.
#' @return \code{list} with indices \code{ik} of extreme values \code{minmax}.
#' @family mbqn
# #' @concept quantile, quantile normalization, rank invariance
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3) # Create a data matrix
#' get_kminmax(X, k = 5, "max") # get indices of the 5 largest values in each column
#' @author Ariane Schad
#  Aug. 2017
#' @export get_kminmax
get_kminmax <- function(X,k,flag = "max"){

  if(flag == "min"){
    decreasing <- FALSE
  }else if(flag=="max"){
    decreasing <- TRUE
  }

  ## sort ##
  na.last <- TRUE
  ik <- apply(X,2,order, decreasing = decreasing, na.last = na.last )
  minmax <- apply(X,2,function(x)x[order(x, decreasing = decreasing, na.last = na.last )])

  ## select the first k elements ##
  ik <- ik[1:k,]
  minmax <- minmax[1:k,]

  ## return results --
  return(list(ik = ik, minmax = minmax))
}
