#' Get the k largest/smallest values and their indices for each
#' column of a matrix
#'
#' @description Extract the k largest/smallest values and their indices for each column of a matrix.
#' @param X data matrix.
#' @param k number of selected extreme values. Must be less then #rows of X.
#' @param flag "min or "max"
#' @details Function used by \code{mbqn.check_saturation}. Sort values of
#' each column of a matrix or an array X and determine the
#' k smallest (flag <- "min") or largest (flag <- "max") values in. NA's in the
#' data are ignored.
#' @return \code{list} with indices \code{ik} of extreme values \code{minmax}.
# #' @keywords quantile normalization, proteomics
#' @family mbqn
#' @concept quantile, quantile normalization, rank invariance
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3) # Create a data matrix
#' get_kminmax(X, k = 5, "max") # get indices of the 5 largest values in each column
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#  Aug. 2017
#' @export get_kminmax
get_kminmax <- function(X,k,flag){

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
