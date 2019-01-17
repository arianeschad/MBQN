#' Selective mean-balanced quantile normalization for (nearly) rank invariant features
#'
#' @description Quantile normalization of a data matrix where selected rows or rank invariant (RI)/
#' nearly rank invariant (NRI) rows/features are normalized by
#' the mean/median balanced quantile normalization.
# Quantile normalization of a data matrix where selected rows or rows with a rank
# invariance frequency above a threshold (=NRI features) are normalized with
# the mean balanced quantile normalization.
#' @inheritParams mbqn
#' @param x A data matrix, where rows represent proteins and
#' columns samples from different replicates, treatments, or conditions.
# #' @param FUN A function like mean, median, a user defined function, or a numeric array
# #' of weights with \code{dim(user_array) = nrow(x)} to balance each intensity profile across samples.
# #' Functions can be parsed also as characters. If FUN = NULL, features are not balanced,
# #' i.e. normal QN is used.
# #' @param method Method used to compute quantile normalization; default NULL - use QN
#' function from the preprocessCore package ; if "limma" - the QN function from the
#' Limma package is used.
# #' @param na.rm A logical value indicating whether NA values should be omitted in the
#' computation of average feature expression.
#' @param index A single index or a vector of indices of selected rows.
#' @param low_thr Value between \[0 1\] as lower threshold that specifies the critical
#' rank invariance frequency, with default set to 0.5.
# #' @param ... Optional arguments passed to \code{mbqn()}
#' @details Specified rows and/or rows with a rank invariance frequency between samples above \code{threshold}
#' are normalized with the mean/median balanced quantile normalization. Remaining rows are quantile normalized
#' without balancing of means.
#' @return Normalized \code{matrix}.
#' @importFrom stats median sd
# #' @keywords Modified Quantile normalization, proteomics.
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples ## Quantile normalize a data matrix where
#' ## nearly rank invariant (NRI) features are balanced
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqn.nri(X, median,low_thr = 0.5) # Balance NRI features selected by threshold
#' mbqn.nri(X, median, index = c(1,2)) # Balance selected features
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export mbqn.nri
# Created: Nov 2018

mbqn.nri <- function(x, FUN = NULL, na.rm = TRUE, method = NULL, low_thr = 0.5 ,index = NULL, verbose = TRUE){

  #opt.args <- list(...)
  #verbose <- ifelse(!is.null(opt.args$verbose), opt.args$verbose, TRUE)

    if(is.null(index)){
        if (!is.numeric(low_thr) || low_thr >1)
          { stop("Wrong data format for low_thr, must be a value within [0 1]!")}

      res  <- mbqn.get_nri_features(x, FUN = FUN,
                                    low_thr = low_thr,
                                    method = method,
                                    verbose = verbose)
      index <- as.numeric(names(res$nri))
      }
  x.mbqn <- mbqn(x = x, FUN = FUN,na.rm = na.rm, method = method, verbose = verbose)
  x.qn <- mbqn(x = x, FUN = NULL ,na.rm = na.rm, method = method, verbose = verbose)
  if(length(index)>0) x.qn[index,] <- x.mbqn[index,]
return(x.qn)
}
