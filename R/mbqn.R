#' Mean/Median-balanced quantile normalization
#'
#' @description Modified quantile-normalization (QN) of a matrix, e.g.,
#' intensity values from omics data or other data sorted in columns.
#' The modification prevents systematic flattening of features (rows) which are
#' rank invariant (RI) or nearly rank invariant (NRI) across columns, for
#' example features that populate mainly the tails of the intensity
#' distribution or features that separate in intensity.
#' @param x a data matrix, where rows represent features, e.g. of protein
#' abundance, and columns represent groups or samples, e.g. replicates,
#' treatments, or conditions.
#' @param FUN a function like mean, median (default), a user defined function,
#' or a numeric vector of weights with length \code{nrow(x)} to balance each
#' feature across samples. Functions can be parsed also as characters.
#' If FUN = NULL, features are not balanced, i.e. normal QN is used.
#' @param na.rm logical indicating to omit NAs in the
#' computation of feature mean.
#' @param method character specifying function for computation of quantile
#' normalization; "limma" (default) for \code{normalizeQuantiles()} from the
#' limma package or "preprocessCore" for \code{normalize.quantiles()} from the
#' preprocessCore package.
#' @param offsetmatrix logical indicating if offset matrix should be used 
#' instead of offset vector specifying offset for each row
#' @param verbose logical indicating to print messages.
#' @details Balance each matrix row by substracting its feature offset computed
#' with FUN, e.g. the median; apply quantile-normalization and add the feature
#' means to the normalized matrix.
#' For further details see \[4\]. For quantile normalization with the "limma"
#' package see \[1,2\] and for the preProcessCore package see \[3\].
#' @return Normalized matrix
#' @importFrom limma normalizeQuantiles
#' @seealso [mbqnNRI()], [mbqnGetNRIfeatures()].
#' @references
#' \[1\] Smyth, G. K., and Speed, T. P. (2003). Normalization of cDNA microarray
#' data. Methods 31, 265–273. \cr
#' Ritchie, M.E., Phipson, B., Wu, D., Hu, Y., Law, C.W., Shi, W., and Smyth,
#' \[2\] G.K. (2015). limma powers differential expression analyses for
#' RNA-sequencing and microarray studies. Nucleic Acids Research 43(7), e47.\cr
#' \[3\] Bolstad B. M. (2016). preprocessCore: A collection of pre-processing
#' functions. R package version 1.36.0.
#' https://github.com/bmbolstad/preprocessCore \cr
#' \[4\] Brombacher, E., Schad, A., Kreutz, C. (2020). Tail-Robust Quantile 
#' Normalization. BioRxiv.
#' @examples
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
#' @author Ariane Schad
#' @export mbqn
# Created: July 2017

mbqn <- function(x, FUN = "mean", na.rm = TRUE, method = "limma", 
    offsetmatrix = FALSE, verbose = FALSE){

    if (is.null(method)) method <- "limma"
    # Check if package limma is installed to run this function
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("Package \"pkg\" needed for this function to work. 
        Please install it.",
        call. = FALSE)
    }

    if (method == "preprocessCore"){
        # Check if package preprocessCore is installed  to run this function
        if (!requireNamespace("preprocessCore", quietly = TRUE)) {
            stop("Package \"pkg\" is required for this function to 
            work. Please install it.",
            call. = FALSE)}
    }

    if (!is.matrix(x)) {
        stop("Wrong data format! Input must be a matrix!")
    }

    # check if data contains NaN and replace it with NA,
    # since preprocessCore will give erronous results in this case
    if (length(which(is.nan(x)))>0)
        x[is.nan(x)] <- NA
    
    if (offsetmatrix) 
        x.napattern <- !is.na(x)

    if (!is.null(FUN)){
        if (is.character(FUN)) FUN <- match.fun(FUN)
        if (is.function(FUN)){
            mx <- apply(x,1,FUN,na.rm=na.rm) # row mean
            
            if (offsetmatrix){
                offset.matrix <-matrix(rep(mx, times=ncol(x)), ncol=ncol(x))
                offset.matrix[x.napattern==FALSE] <- NA
                
                if (is.null(method) || method == "limma"){
                    mx <- normalizeQuantiles(offset.matrix)
                    rownames(mx) <- NULL
                } else if (method == "preprocessCore"){
                    mx <- preprocessCore::normalize.quantiles(offset.matrix)
                }
            }
        }
        if (is.numeric(FUN)){
            mx <- FUN
            if (sum(abs(FUN)==0, na.rm =TRUE)){
                message("Array-elements are all zero. Compute QN 
                without mean balancing.")
            }
        }

        # balanced quantile normalisation
        if (is.null(method) || method == "limma"){
            dummy <- normalizeQuantiles(x-mx)
            rownames(dummy) <- NULL
        } else if (method == "preprocessCore"){
            dummy <- preprocessCore::normalize.quantiles(x-mx)
        }
        qn_x <- dummy + mx

    } else {
        if (verbose) message("Compute QN without mean balancing.")
        # quantile normalisation
        if (is.null(method) || method == "limma") {
            qn_x <- normalizeQuantiles(x)
            rownames(qn_x) <- NULL
        } else if (method == "preprocessCore"){
            qn_x <- preprocessCore::normalize.quantiles(x)
        }
    }
    colnames(qn_x) <- colnames(x)
    row.names(qn_x) <- row.names(x)
    return(qn_x)
}
