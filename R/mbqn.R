#' Mean-balanced quantile normalization
#'
#' @param x A matrix where rows represent features, e.g. protein abundances/intensities and
#' columns are samples from replicates or conditions.
#' @param FUN A function like mean, median, or a user defined function or an array with
#' \code{dim(user_array) = nrow(x)}. Default NULL - features are not balanced.
#' @param method Function to compute quantile normalization; default NULL - use function from the preprocessCore package ; if "limma" - the function from the Limma package is used.
#' @export
#' @details Normalize a data matrix based on a mean-balanced quantile normalization.
#' Each row of the data matrix is balanced by FUN, e.g. the median, before normalization.
#' After normalization, row means are added to the normalized matrix.
#' This function uses \code{preprocessCore::normalize.quantiles()} by Bolstad et al, Bioinformatics (2003),
#' available from http://bitools::package_dependencies(pkgs, db, which = 'all', reverse = TRUE)bioconductor.org/biocLite.R.
#  see source('http://bioconductor.org/biocLite.R') and biocLite('preprocessCore').
#' @return Mean-/Median-balanced quantile normalized \code{matrix}.
#' @keywords Modified Quantile normalization,  proteomics.
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
# Schad, A. and Kreuz, C. (2017) Mean-balanced quantile
# normalization for processing label-free quantitative proteomics
# data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples ## Compute quantile normalization using preprocessCore
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3)
#' mbqn(X)
#'
#' ## Compute median and mean balanced quantile normalization
#' mbqn(X, median) # Use median to center features
#' mbqn(X, mean) # Use mean to center features
#'
#' ## Use user defined array of weighted averages for centering
#' wt <- c(1,3,1)/5 # Weights for each sample
#' user_array <- apply(X,1,weighted.mean, wt ,na.rm =T)
#' mbqn(X, user_array)
#'
#'## Use limma package to compute quantile normalization
#' mbqn(X, median, method = "limma")
#' @description Modified quantile-normalization of a matrix, representing for example
#' omics or other data sorted in a matrix. Prevents systematic flattening of feature variation across columns
#' for features overrepresented in the tails of the intensity distribution
#' across columns, i.e. rank invariant (RI) or nearly rank invariant (NRI) features.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export
# Created: July 2017 - update Oct. 2018

mbqn <- function(x, FUN = NULL, na.rm = TRUE, method = NULL){

  # Check if package preprocessCore is installed  to run this function
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package \"pkg\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

  if(is.character(FUN)) FUN <- match.fun(FUN)

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
    if(!is.null(method) && method == "limma"){
      dummy <- limma::normalizeBetweenArrays(x-mx)
      rownames(dummy) <- NULL
    }else{
      dummy <- preprocessCore::normalize.quantiles(x-mx)
      }
      qn_x <- dummy + mx

  } else {
    print("Comput QN without mean balancing.")
    # quantile normalisation
    if(!is.null(method) && method == "limma") {
      qn_x <- limma::normalizeBetweenArrays(x)
      rownames(dummy) <- NULL
    }else{
      qn_x <- preprocessCore::normalize.quantiles(x)
    }
  }
  return(qn_x)
}
