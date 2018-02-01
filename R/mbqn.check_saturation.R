#' Test if quantile-normalization of data benefits from mean-balancing
#'
#' Rank data and test if lower and upper intensity tails are dominated by few feature types/proteins. Compute a quantile normalization without and with mean-balancing and check standard deviation of normalized data entries located in the tails
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param mean_fun mean or median, if left empty, quantile normalization is applied without balancing the data
#' @param qlow lower quantile
#' @param qup upper quantile, default 1
#' @return A matrix of mean- or median-balanced quantile normalized data
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile normalization for processing label-free quantitative proteomics data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples mbqn(dat, mean)
#' @description This function uses normalize.quantiles() from the package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
#' @author A. Schad, Aug. 2017


mbqn.check_saturation <- function(dat, FUN = mean_fun, qlow, qup){

  #quantile normalisation
  qn.dat <- preprocessCore::normalize.quantiles(dat)
  mbqn(dat,FUN = mean, na.rm = TRUE)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)
  # sample mean for each row (protein)
  m.dat <- apply(dat,1,FUN, na.rm=TRUE)
  #mean balanced quantile normalisation
  mbqn.dat <- preprocessCore::normalize.quantiles(dat-mdat) + mdat
  s.mbqn <- apply(mbqn.dat, 1, sd, na.rm=TRUE)

  return(qn_dat)
}
