#' Mean-balanced quantile normalization
#'
#' @param x a data matrix, where rows represent intensity features, e.g.
#' of protein abundance, and columns represent samples, e.g. replicates, treatments, or conditions.
#' @param FUN a function like mean, median, a user defined function, or a numeric array
#' of weights with
#' \code{dim(user_array) = nrow(x)} to balance each intensity profile across samples.
#' Functions can be parsed also as characters. If FUN = NULL, features are not balanced,
#' i.e. normal QN is used.
#' @param na.rm a logical value indicating whether NA values should be omitted in the
#' computation of average feature expression.
#' @param method packagename providing the function used to compute the standard quantile normalization;
#' default NULL or "limma" - use \code{normalizeBetweenArrays()} from the limma package,
#' if "preprocessCore" - \code{normalize.quantiles()} from the preprocessCore package is used.
#' @param verbose a logical indicating to run function quietly.
#' @details Normalize a data matrix based on a mean-balanced quantile normalization.
#' Each row of the data matrix is balanced by FUN, e.g. the median, before normalization.
#' After normalization, row means are added to the normalized matrix.
#' By default, this function uses \code{limma::normalizeQuantiles()}.
# # \code{limma::normalizeBetweenArrays()}
#' available from http://bioconductor.org/biocLite.R.
#' @return Normalized matrix
#' @importFrom limma normalizeQuantiles
# #' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep., 2019
#' @examples
#'\dontrun{
#' ## Compute mean and median balanced quantile normalization
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqn(X, mean) # Use arithmetic mean to center features
#' mbqn(X, median) # Use median to center features
#' mbqn(X, "median")
#'
#' ## Use user defined array of weights for averaging
#' wt <- c(1,3,1)/5 # Weights for each sample
#' user_array <- apply(X,1,weighted.mean, wt ,na.rm =TRUE)
#' mbqn(X, user_array)
# #'
# #' ## Use limma package to compute quantile normalization
# #' mbqn(X, median, method = "preprocessCore")
#' }
#' @description Modified quantile-normalization of a matrix, representing for example
#' omics or other data sorted in a matrix. Prevents systematic flattening of feature variation across columns
#' for features overrepresented in the tails of the intensity distribution
#' across columns, i.e. rank invariant (RI) or nearly rank invariant (NRI) features.
#' @author Ariane Schad
#' @export mbqn
# Created: July 2017

mbqn <- function(x, FUN = NULL, na.rm = TRUE, method = NULL, verbose = TRUE){

  # Check if package limma is installed to run this function
  if (!requireNamespace("limma", quietly = TRUE)) {
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
    if(is.null(method) || method == "limma"){
      dummy <- normalizeQuantiles(x-mx)
      rownames(dummy) <- NULL
    } else if (method == "preprocessCore"){
      # Check if package preprocessCore is installed  to run this function
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        stop("Package \"pkg\" needed for this function to work. Please install it.",
             call. = FALSE)}
      dummy <- preprocessCore::normalize.quantiles(x-mx)
    }
    qn_x <- dummy + mx

  }else{
    if(verbose) print("Comput QN without mean balancing.")
    # quantile normalisation
    if(is.null(method) || method == "limma") {
      qn_x <- normalizeQuantiles(x)
      rownames(qn_x) <- NULL
    }else if (method == "preprocessCore"){
      # Check if package preprocessCore is installed  to run this function
      if (!requireNamespace("preprocessCore", quietly = TRUE)) {
        stop("Package \"pkg\" needed for this function to work. Please install it.",
             call. = FALSE)}
      qn_x <- preprocessCore::normalize.quantiles(x)
    }
  }
  colnames(qn_x) <- colnames(x)
  return(qn_x)
}
