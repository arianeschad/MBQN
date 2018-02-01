#' Mean-balanced quantile normalization
#'
#' Function to compute a mean-balanced quantile normalization.
#' @param x A data matrix. Rows represent features, e.g. protein abundances; columns replicates, samples.
#' @param FUN mean, median, or a user defined function, or an array with dim = numbers of rows of dat. If left empty, quantile normalization is applied without balancing the data.
#' @return A matrix of mean- or median-balanced quantile normalized data.
#' @keywords quantile normalization proteomics.
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile normalization for processing label-free quantitative proteomics data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples mbqn(x, mean)
#' mbqn(x, median)
#' mbqn(x, user_function)
#' mbqn(x, user_array)
#' @description This function uses normalize.quantiles() from the package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
#' @author A. Schad, Aug. 2017


mbqn <- function(x, FUN = NULL, na.rm = TRUE){

  if (!is.matrix(x)) {
    stop("Wrong data format! Input x must be a matrix!")
  }

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
