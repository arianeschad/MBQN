#' Identify rank invariant (RI) and nearly rank invariant (NRI) features
#'
#' @param dat a data matrix. Rows represent features, e.g. protein abundances; columns represent samples
# #' @param FUN median, mean, or another function used to balance features across colums. If left empty, quantile normalization
#' is applied without balancing the data
#' @inheritParams mbqnNRI
#' @importFrom stats median sd
#' @description Compute the rank frequency of each feature of a matrix and identify NRI/RI features.
# #' @return A list of identified RI and NRI features, rank frequencies, and rank variation. p = p, max_p = max_p, ip = ip, nri = nri, var0_feature = ind_var0
#' @return A list with elements:
#' \item{\code{p}}{a matrix of the rank invariance frequencies and the sample coverage over all RI/NRI features}
#' \item{\code{max_p}}{value of the maximum rank invariance frequency in percent}
#' \item{\code{ip}}{index of the feature with maximum rank invariance frequency}
#' \item{\code{nri}}{a table of the rank invariance frequencies in percent of NRI/RI features}
#' \item{\code{var0_feature}}{index of features with zero sample variance after QN.}
# #' @keywords quantile normalization, proteomics
# #' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
# @examples ## Check data matrix for RI and NRI features
# X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
# mbqn.check_saturation(X, mean, low_thr = 0.5, save_fig = FALSE)
#' @details Quantile normalize a data matrix, sort ranks, and count the maximum times features share
#' the same ranks across all columns. A user-defined threhold is used to define nearly rank invariant features.
#' @author Ariane Schad
#' @export mbqnGetNRIfeatures
#  Created: Nov 2018
mbqnGetNRIfeatures <- function(dat, FUN = NULL,
                               low_thr = 0.5,
                               method = NULL,
                               verbose = TRUE){

  # if FUN is not specified, use median!
  if(is.null(FUN)) FUN <- median
  if(is.character(FUN)) FUN <- match.fun(FUN)

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL, method = method, verbose = verbose)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  ## Rank frequencies for each feature after QN (top-down)
  # & assign NAs to 0 rank

  out <- MBQN::get_kminmax(X = qn.dat,k = N, flag = "max")

  tdummy = lapply(
    1:N,
    function(i)
    {
      table(which(out$ik==i,arr.ind =TRUE)[,1]*(!is.na(dat[i,])),useNA = 'ifany')
    }
  )

  # scale to frequencies
  pi <- lapply(tdummy,function(x) x/M)

  max_pi = lapply(
    1:N,
    function(i)
    {
      # ignore NAs!
      ki <- bla <- NULL
      ki <- which(as.numeric(names(pi[[i]]))!=0)
      if(length(ki)>0){
        bla <- pi[[i]][which(pi[[i]]==max(pi[[i]][ki], na.rm =TRUE))]
      }else{
        bla <- pi[[i]][1]*0
        #bla[[1]][1] <- 1 # if it has no rank
      }
    }
  )

  max_pi_vals <- lapply(max_pi,function(i) unique(i))
  max_pi_vals[which(sapply(max_pi_vals,length)<1)] <- 0

  # how often are data present for these features
  not_nas <- apply(qn.dat,1,function(x) length(which(!is.na(x))))/M

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

    freq_ismissing = sum(is.na(qn.dat[ip,]))/dim(qn.dat)[2] # cnt how often protein is missing

    if(verbose) print(paste('Maximum frequency of RI/NRI feature(s): ',max_p,"%"))

    #########################################################################################

    nri <- p*100

    p <- rbind(p,not_nas[as.numeric(names(p))])
    rownames(p) <- c("occupation.freq", "sample.coverage")

    # which features have zero variation after QN
    ind_var0 <- which(s.qn==0)
    } else {
    if(verbose) print(paste('No RI/NRI feature(s) found!'))
    max_p <- NULL
    ip <- NULL
    nri <- NULL

  }
  ind_var0 <- which(s.qn==0)

  return(list(p = p, max_p = max_p, ip = ip, nri = nri, var0_feature = ind_var0))
}
