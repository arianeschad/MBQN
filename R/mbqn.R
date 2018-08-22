#' Mean-balanced quantile normalization
#'
#' @param x A matrix where rows represent features, e.g. protein abundances/intensities and
#' columns are replicates or probes.
#' @param FUN A function like mean or median or a user defined function or an array with
#' \code{dim(user_array) = nrow(x)}. If left empty, the matrix is not balanced.
#' #@export
#' @details Function to normalize a data matrix based on a mean-balanced quantile normalization.
#' Each row of the data matrix is balanced by its mean/median before normalization.
#' Row means are added to the normalized matrix.
#' This function uses \code{preprocessCore::normalize.quantiles()} by Bolstad et al, Bioinformatics (2003),
#' installed from http://bioconductor.org/biocLite.R.
#  by source('http://bioconductor.org/biocLite.R') biocLite('preprocessCore').
#' @return Mean-/Median-balanced quantile normalized \code{matrix}.
#' @keywords Modified Quantile normalization proteomics.
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile
#' normalization for processing label-free quantitative proteomics
#' data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples mbqn(x, mean)
#' mbqn(x, median)
#' mbqn(x, user_function)
#' mbqn(x, user_array)
#' @description Modified quantile-normalization of a matrix representing
#' omics or microarray data. Suppress systematic flattening of feature variation across columns
#' for features overrepresented in the tails of the intensity distribution
#' across columns.
# @seealso \code{\link{xxx}}
#' @author A. Schad, \email{ariane.schad@zbsa.de}
# Aug. 2017
#' @export
mbqn <- function(x, FUN = NULL, na.rm = TRUE){

  if (!is.matrix(x)) {
    stop("Wrong data format! Input x must be a matrix!")
  }

  # check if data contains NaN and replace it with NA, since preprocessCore will
  # give erronous results in this case
  if (length(which(is.nan(x)))>0)
    x[is.nan(x)] <- NA

  if(!is.null(FUN)){
    if(is.function(FUN)){
      mx <- apply(x,1,FUN,na.rm=na.rm) # row mean
    } else if(is.array(FUN)){
      mx <- FUN
    }
    # balanced quantile normalisation
    qn_x <- preprocessCore::normalize.quantiles(x-mx) + mx
  } else {
    print("Comput QN without mean balancing.")
    #quantile normalisation
    qn_x <- preprocessCore::normalize.quantiles(x)
  }
  return(qn_x)
}
