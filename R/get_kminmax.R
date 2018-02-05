#' Get the k largest/smallest values and their indices for each
#' column of a matrix
#'
#' @description Generate a data matrix for illustration of mean-balanced.
#' quantile normalization.
#' @param X data matrix.
#' @param k number of searched extreme values.
#' @param flag "min or "max"
#' @details Subfunction used by \code{mbqn.check_saturation}. Search for
#' k smallest (flag <- "min") or largest (flag <- "max") values in
#' each column of an array or a matrix X. NA's are ignored in the data.
#' @return \code{list} with indices \code{ik} of extreme values \code{minmax}.
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced
#' quantile normalization for processing label-free quantitative
#' proteomics data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples
#' get_kminmax(X, k = 5, "max")
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#
# A. Schad, Aug. 2017

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

  ## return results
  return(list(ik = ik, minmax = minmax))
}
