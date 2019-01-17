#' A core function of MBQN. Identify and extract frequencies of rank invariant (RI) and nearly rank invariant (NRI) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
# #' @param FUN median, mean, or another function used to balance features across colums. If left empty, quantile normalization
#' is applied without balancing the data
# #' @param qlow lower quantile
# #' @param qup upper quantile, default 1
# #' @param low_thr Numerical value for the lower threshold for NRI frequency, default = 0.2
# #' @param method Packagename containing function used to compute quantile normalization; default NULL - use the preprocessCore package ; "limma" uses the Limma package
# #' @param verbose Logical for running function quiet
#' @inheritParams mbqn.nri
#' @details Compute rank frequencies rank variation of each feature of a quantile normalized matrix.
#' @return List of identified RI and NRI features, rank frequencies, and rank variation.
# #' @keywords quantile normalization, proteomics
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
# @examples ## Check data matrix for RI and NRI features
# X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
# mbqn.check_saturation(X, mean, low_thr = 0.5, save_fig = FALSE)
#' @description Apply quantile normalization to data matrix, sort ranks and count times features share
#' same ranks across samples/columns. A user-defined threhold is used to identify nearly rank invariant features.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export mbqn.get_nri_features
#  Created: Nov 2018
mbqn.get_nri_features <- function(dat, FUN = NULL,
                                  low_thr = 0.5,
                                  method = NULL,
                                  verbose = TRUE){

  # if FUN is not specified, use median!
  if(is.null(FUN)) FUN <- median
  if(is.character(FUN)) FUN <- match.fun(FUN)

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL, verbose = verbose)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  # quantile normalisation and its standard deviation computed with the limma function
  # The results from limma::normalizeBetweenArrays normalization are equal, up to numerical differences,
  # to that of preprocessCore::normalize.quantiles!

  ##############################################################################
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
      ki <- which(as.numeric(names(pi[[i]]))!=0)
      bla <- pi[[i]][which(pi[[i]]==max(pi[[i]][ki], na.rm =TRUE))]
    }
  )

  max_pi_vals <- lapply(max_pi,function(i) unique(i))

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
    # how often are data present for these features
    # if(length(ind_var0)>0){
    #    not_nas <- apply(qn.dat,1,function(x) length(which(!is.na(x))))/M
    # not_nas <- not_nas[ind_var0]
    # names(not_nas) <- ind_var0
    # not_nas <- as.table(not_nas)
    # } else {not_nas <- NULL}
  } else {
    if(verbose) print(paste('No RI/NRI feature(s) found!'))
    max_p <- NULL
    ip <- NULL
    nri <- NULL

  }
  ind_var0 <- which(s.qn==0)

  return(list(p = p, max_p = max_p, ip = ip, nri = nri, var0_feature = ind_var0))
}
