#' Selective mean-balanced quantile normalization for (nearly) rank invariant features
#'
#' @description Quantile normalization of a data matrix where selected rows or rank invariant (RI)/
#' nearly rank invariant (NRI) rows/features are normalized by
#' the mean/median balanced quantile normalization.
# Quantile normalization of a data matrix where selected rows or rows with a rank
# invariance frequency above a threshold (=NRI features) are normalized with
# the mean balanced quantile normalization.
#' @inheritParams mbqn
#' @param index A single index or a vector of indices of selected rows.
#' @param low_thr Value between \[0 1\] as lower threshold that specifies the critical
#' rank invariance frequency, with default set to 0.5.
#' @details Specified rows and/or rows with a rank invariance frequency between samples above \code{threshold}
#' are normalized with the mean/median balanced quantile normalization. Remaining rows are quantile normalized
#' without balancing of means.
#' @return Normalized \code{matrix}.
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples ## Quantile normalize a data matrix where
#' ## nearly rank invariant (NRI) features are balanced
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqnNRI(X, median,low_thr = 0.5) # Balance NRI features selected by threshold
#' mbqnNRI(X, median, index = c(1,2)) # Balance selected features
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export mbqnNRI
# Created: Nov 2018

mbqnNRI <- function(x, FUN = NULL, na.rm = TRUE, method = NULL, low_thr = 0.5 ,index = NULL, verbose = TRUE){

  if(is.null(index)){
    if (!is.numeric(low_thr) || low_thr >1)
    { stop("Wrong data format, low_thr must be a value within [0 1]!")}

    res  <- mbqnGetNRIfeatures(x, FUN = FUN,
                               method = method,
                               low_thr = low_thr,
                               verbose = verbose)
    index <- as.numeric(names(res$nri))
  }
  x.mbqn <- mbqn(x = x, FUN = FUN,na.rm = na.rm, method = method, verbose = verbose)
  x.qn <- mbqn(x = x, FUN = NULL ,na.rm = na.rm, method = method, verbose = verbose)
  if(length(index)>0) x.qn[index,] <- x.mbqn[index,]
  return(x.qn)
}
