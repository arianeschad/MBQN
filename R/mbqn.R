#' Mean-balanced quantile normalization
#'
#' @param x A data matrix, where rows represent proteins and
#' columns samples from different replicates, treatments, or conditions.
#' @param FUN A function like mean, median, or a user defined function, or a numeric array with
#' \code{dim(user_array) = nrow(x)} to balance each intensity profile across samples.
#' Functions can be parsed also as characters. If FUN = NULL, features are not balanced, i.e. normal QN is used.
#' @param method Function used to compute quantile normalization; default NULL - use function from the preprocessCore package ; if "limma" - the function from the Limma package is used.
#' @param na.rm A logical value indicating whether NA values should be omitted in the computation of average feature expression.
#' @param verbose Logical to print annotations to console
#' @details Normalize a data matrix based on a mean-balanced quantile normalization.
#' Each row of the data matrix is balanced by FUN, e.g. the median, before normalization.
#' After normalization, row means are added to the normalized matrix.
#' This function uses \code{preprocessCore::normalize.quantiles()} by Bolstad et al, Bioinformatics (2003),
# #' available from http://bitools::package_dependencies(pkgs, db, which = 'all', reverse = TRUE)bioconductor.org/biocLite.R.
#' available from http://bioconductor.org/biocLite.R.
# # from source('http://bioconductor.org/biocLite.R') and biocLite('preprocessCore').
#' @return Mean-/Median-balanced quantile normalized \code{matrix}.
#' @importFrom stats median sd
# #' @keywords Modified Quantile normalization, proteomics.
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note, 2018
# Schad, A. and Kreuz, C. (2017) Mean-balanced quantile
# normalization for processing label-free quantitative proteomics
# data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples ## Compute quantile normalization using preprocessCore
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqn(X)
#'
#' ## Compute mean balanced quantile normalization
#' mbqn(X, mean) # Use arithmetic mean to center features
#' mbqn(X, median) # Use median to center features
#' mbqn(X, "median") # Use median to center features
#'
#' ## Use user defined array of weights for averaging
#' #wt <- c(1,3,1)/5 # Weights for each sample
#' #user_array <- apply(X,1,weighted.mean, wt ,na.rm =TRUE)
#' #mbqn(X, user_array)
#'
#'## Use limma package to compute quantile normalization
#' mbqn(X, median, method = "limma")
#' @description Modified quantile-normalization of a matrix, representing for example
#' omics or other data sorted in a matrix. Prevents systematic flattening of feature variation across columns
#' for features overrepresented in the tails of the intensity distribution
#' across columns, i.e. rank invariant (RI) or nearly rank invariant (NRI) features.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export mbqn
# Created: July 2017 - update Oct. 2018

mbqn <- function(x, FUN = NULL, na.rm = TRUE, method = NULL, verbose = TRUE){

  # Check if package preprocessCore is installed  to run this function
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package \"pkg\" needed for this function to work. Please install it.",
           call. = FALSE)
    }

  if (!is.matrix(x)) {
    stop("Wrong data format! Input must be a matrix!")
  }

  # check if data contains NaN and replace it with NA, since preprocessCore will
  # give erronous results in this case
  if (length(which(is.nan(x)))>0)
    x[is.nan(x)] <- NA

  if(!is.null(FUN)){
    if(is.character(FUN)) FUN <- match.fun(FUN)
    if(is.function(FUN)){
      mx <- apply(x,1,FUN,na.rm=na.rm) # row mean
      }
    if(is.numeric(FUN)){mx <- FUN
    if(sum(abs(FUN)==0,na.rm =T))  print("Array-elements are all zero. Comput QN without mean balancing.")
    }
    # balanced quantile normalisation
    if(!is.null(method) && method == "limma"){
      dummy <- limma::normalizeBetweenArrays(x-mx)
      rownames(dummy) <- NULL
    }else{
      dummy <- preprocessCore::normalize.quantiles(x-mx)
      }
      qn_x <- dummy + mx

  }else{
    if(verbose) print("Comput QN without mean balancing.")
    # quantile normalisation
    if(!is.null(method) && method == "limma") {
      qn_x <- limma::normalizeBetweenArrays(x)
      rownames(dummy) <- NULL
    }else{
      qn_x <- preprocessCore::normalize.quantiles(x)
    }
  }
  colnames(qn_x) <- colnames(x)
  return(qn_x)
}
