#' Identify rank invariant (RI) and nearly rank invariant (NRI) features
#'
#' @description Compute the rank frequency of each feature of a matrix and
#' identify NRI/RI features.
#' @param x a data matrix. Rows represent features, e.g. protein abundances;
#' columns represent samples.
#' @param low_thr a value between \[0 1\]. Features with RI
#' frequency >=\code{low_thr} are considered as NRI/RI; default 0.5.
#' @inheritParams mbqn
#' @importFrom stats median sd
#' @return A list with elements:
#' \item{\code{p}}{a matrix with the rank invariance frequencies \code{ri.freq}
#' and the sample coverage \code{sample.coverage} for all detected RI/NRI
#' features}
#' \item{\code{max_p}}{maximum rank invariance frequency in percent}
#' \item{\code{ip}}{index of feature with maximum rank invariance frequency}
#' \item{\code{nri}}{table of the rank invariance frequencies in percent for
#' each NRI/RI feature}
#' \item{\code{var0_feature}}{indices of features with zero sample variance
#' after QN}
#' \item{\code{low_thr}}{threshold used for NRI/RI detection from RI frequency.}
#' @seealso [mbqnPlotRI()] for visualization of detected NRI/RI features.
#' @references Brombacher, E., Schad, A., Kreutz, C. (2020). Tail-Robust 
#' Quantile Normalization. BioRxiv.
#' @examples ## Check data matrix for RI and NRI features
#' set.seed(1234)
#' x <- mbqnSimuData("omics.dep")
#' RI <- mbqnGetNRIfeatures(x, low_thr = 0.5, verbose = FALSE)
#' mbqnPlotRI(RI)
#' @details Quantile normalize the data matrix and sort ranks. Determine the
#' maximum frequency of equal rank across all columns for each feature.
#' Features with maximum frequency above the user-defined threhold are declared
#' as nearly rank invariant.
#' @author Ariane Schad
#' @export mbqnGetNRIfeatures
#  Created: Nov 2018
mbqnGetNRIfeatures <- function(x,
  low_thr = 0.5,
  method = NULL,
  verbose = TRUE){
  
  N <- nrow(x)
  M <- ncol(x)
  
  # classical quantile normalisation and its standard deviation
  qn.x <- mbqn(x = x, FUN = NULL, method = method, verbose = FALSE)
  s.qn <- apply(qn.x, 1, sd, na.rm=TRUE)
  
  ## Rank frequencies for each feature after QN (top-down)
  # & assign NAs to 0 rank
  
  out <- getKminmax(x = qn.x, k = N, flag = "max")
  
  tdummy = lapply(
    seq_len(N),
    function(i){
      # set rank to zero for NA and count the number of each rank per row
      table(which(out$ik==i, arr.ind =TRUE)[,1]*(!is.na(x[i,])),
        useNA = 'ifany')
    })
  
  range01 <- function(z){(z-min(z))/(max(z)-min(z))}
  
  tdummy2 <- lapply(tdummy, function(y) {
    ki <- which(as.numeric(names(y))==0)
    if (length(ki) > 0){
      y <- y[-ki]
    }
    y
  })
  
  # Scale to frequencies (without taking NAs into account)
  pi <- lapply(tdummy2,function(x) x/sum(x))
  
  max_pi = lapply(
    seq_len(N),
    function(i){
      # Ignore NAs!
      ki <- maxfreq <- NULL
      ki <- which(as.numeric(names(pi[[i]]))!=0)
      if(length(ki)>0){
        maxfreq <- pi[[i]][which(pi[[i]]==max(pi[[i]][ki], na.rm =TRUE))]
      }else{
        maxfreq <- pi[[i]][1]*0
      }
    })
  
  max_pi_vals <- lapply(max_pi,function(i) unique(i))
  
  max_pi_vals[which(vapply(max_pi_vals,length,FUN.VALUE = numeric(1))<1)] <- 0
  
  # how often are data present for these features
  not_nas <- apply(qn.x,1,function(x) length(which(!is.na(x))))/M
  
  p <- max_pi_vals[which(max_pi_vals>=low_thr)]
  p <- unlist(p)
  
  is.defined <- (!is.null(p))
  
  if(is.defined){
    names(p) <- as.character(which(max_pi_vals>=low_thr))
    p <- as.table(p)
    
    if(length(p)>1){
      if(verbose) message('Caution: There might be multiple RI/NRI features!')
    }
    
    max_p <- max(p)*100
    ip <- as.integer(names(which(p*100==max_p)))
    
    # Count how often protein is missing
    freq_ismissing = sum(is.na(qn.x[ip,]))/ncol(qn.x)
    
    if (verbose) message(paste('Maximum frequency of RI/NRI feature(s): ',
      max_p,"%"))
    
    # In percent
    nri <- p*100
    
    p <- rbind(p,not_nas[as.numeric(names(p))])
    rownames(p) <- c("RI.freq", "sample.coverage")
    
    # Which features have zero variation after QN
    ind_var0 <- which(s.qn==0)
  } else {
    if(verbose) message(paste('No RI/NRI feature(s) found!'))
    max_p <- NULL
    ip <- NULL
    nri <- NULL
  }
  ind_var0 <- which(s.qn==0)
  
  return(list(p = p, max_p = max_p, ip = ip, nri = nri,
    var0_feature = ind_var0, low_thr = low_thr))
}