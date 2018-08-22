#' Check data matrix for rank invariant (ri) /nearly rank invariant (nri) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param mean_fun median or mean, if left empty, quantile normalization
#' is applied without balancing the data
#' @param qlow lower quantile
#' @param qup upper quantile, default 1
#' @param flag_show_fig flag to specify whether results are plotted to figure, default = TRUE
#' @inheritParams mbqn
#' @details Rank data and check if lower and upper intensity tails are
#' dominated by few feature. Compute a quantile
#' normalization without and with mean-balancing and check standard
#' deviation of normalized data entries located in the tails. For each feature, count maximum
#' rank frequencies if feature is not dominated by missing values
#' @return A matrix of median- or mean-balanced quantile normalized data
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Median-balanced quantile
#' normalization for processing label-free quantitative proteomics
#' data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples mbqn.check_saturation(dat, mean)
#' @description Test if few rows of the data matrix dominate the upper tail. This script uses normalize.quantiles() from the
#' package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' 2017
#' @export
mbqn.check_saturation_v2 <- function(dat, FUN = mean_fun, qlow, qup, flag_show_fig = TRUE){

  # dat <- matrix(rnorm(20000,0,1),nrow = 2000, ncol = 10)
  # dat[1,] <- dat[1,]+4
  # dat[2,] <- dat[2,]+2
  # dat[3,] <- dat[3,]+1

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # row means
  #m.dat <- apply(dat,1,FUN, na.rm=TRUE)

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  #mean balanced quantile normalisation and its standard deviation
  mbqn.dat <- MBQN:::mbqn(dat,FUN = FUN, na.rm = TRUE)
  s.mbqn <- apply(mbqn.dat, 1, sd, na.rm=TRUE)

#########################################################################################
  # compute (tied) ranks
  ir <- sapply(1:M,function(i) rank(dat[,i], na.last = TRUE, ties.method = c("average")))
  s.ir <- apply(ir, 1, sd, na.rm=TRUE)
  which(s.ir == 0)
  plot(s.ir,type = "l")
  p <- 5

  # find indices of features occupying the top n ranks
  top <- sapply(1:M,function(i) which(ir[,i] %in% c((max(ir)-p):max(ir))))

  # Statistics: shuffle ranks of each column and compute statistics
  s <- matrix(NA,nrow = 20,nco=dim(ir)[1])
  for (i in 1:20){
    test <- sapply(1:dim(ir)[2],function(i) sample(ir[,i],dim(ir)[1],replace = FALSE))
    s[i,] <- apply(test, 1, sd, na.rm=TRUE)
  }
  s.m <- colMeans(s)
  s.s <- apply(s,2,sd,na.rm = TRUE)

  plot(s.ir,type="l",col=c(2))
  lines(s.m)

###########################################################################################
  ## get rank frequencies for each feature after QN (top-down)
  #out <- MBQN:::get_kminmax(X = dat,k = dim(dat)[1], flag = "max")
  #out$ik
  #out$minmax
  #t <- table(out$ik)

  # quantile normalisation
  qn.dat <- MBQN:::mbqn(x = dat,FUN = NULL)

  # mean balanced quantile normalisation
  mbqn.dat <- MBQN:::mbqn(x = dat,FUN = FUN, na.rm = TRUE)

  out=MBQN:::get_kminmax(X = qn.dat,k = 5,flag = "max")
  t = table(out$ik)
  p <- t/M # frequency of feature indices in one of the top k ranks
  max_p <- max(p)
  ip <- as.integer(names(which(t==max(t)))) # index of most frequent top-k level feature

  if(length(ip)>1){
    #ip = ip[1]
    warning('There are multiple RI/NRI features!')
    flag_ambigious<-TRUE
  }
  freq_ismissing = sum(is.na(qn.dat[ip,]))/dim(qn.dat)[2] # cnt how often protein is missing

  print(paste('Maximum frequency of RI/NRI feature(s) in top-k ranks: ',max_p*100,"%"))


  ##################################################################################
  # which features show no variation after QN
  ind_var0 <- which(apply(qn.dat,1,sd, na.rm=T)==0)

  ##############################################################################
  ## get rank frequencies for each feature after QN (top-down)
  # assign NAs to 0 rank

  out <- MBQN:::get_kminmax(X = qn.dat,k = dim(dat)[1], flag = "max")
  wi <- array(NA,dim = dim(qn.dat)[1])
  sw_rel <- max_pi <- wi

  #tdummy <- lapply(1:N,function(i) table(which(out$ik==i,arr.ind =TRUE)[,1],useNA = 'ifany'))
  tdummy <- lapply(1:N, function(i) table(which(out$ik==i,arr.ind =T)[,1]*(!is.na(dat[i,])),useNA = 'ifany'))


  pi <- lapply(tdummy,function(x) x/dim(qn.dat)[2])

  # relative number of ranks populated across samples of each feature
  # small numbers indicated few variation
  nranks <- lapply(pi,function(x) length(x)/M)


  # plot(log(unlist(nranks)))

  #max_pi <- lapply(1:N,function(i) pi[[i]][which(pi[[i]]==max(pi[[i]], na.rm =T))])
  # neglect NAs!
  ki <- which(as.numeric(names(pi[[i]]))!=0)

  max_pi <- lapply(1:N,function(i) pi[[i]][which(pi[[i]]==max(pi[[i]], na.rm =T))])

  max_pi <- lapply(max_pi,function(i) unique(i))


  ind <- which(unlist(max_pi)*100>30)
  pi[ind]

  #bla<- unlist(max_pi[ind])
  #names(bla) <- ind

  # ranks occupied by feature i
  ris <- lapply(tdummy, function(x) as.numeric(names(x)))
  mw <- lapply(1:N,function(i) sum(pi[[i]]*ris[[i]])/sum(pi[[i]]))

  # frequency weighted rank variability around rank frequency weighted mean
  sw_rel <- lapply(1:N,function(i) sqrt(sum(pi[[i]]*(ris[[i]]-mw[[i]])^2)/sum(pi[[i]])))
  max_pi <- unlist(max_pi)
  #plot(unlist(max_pi)*100)
  #ind <- which(max_pi>=0.4)

  #ind <- lapply(1:N, function(i) max_pi[[i]][which(max_pi[[i]]>0.3)])
  ind <- lapply(1:N, function(i) which(max_pi[[i]]>=0.4))
  ind <- unlist(ind)

  #plot(unlist(mw[ind]),max_pi[ind]*100,ylim = c(0,101))
  ########################################################

  # this loop is pretty slow!!!

  # for (i in 1:dim(qn.dat)[1]){
  #   val <- which(out$ik==i, arr.ind = TRUE)
  #   tdummy <- table(val[,1])
  #
  #   # rank frequencies/fraction pi of feature i
  #   pi <- tdummy/dim(qn.dat)[2]
  #   # where(pi,max)
  #   max_pi[i] <- max(tdummy)/dim(qn.dat)[2]
  #   # ranks occupied by feature i
  #   ris <- as.numeric(names(tdummy))
  #   # mean weighted rank
  #   mw[i] <- sum(pi*ris)/sum(pi)
  #   # frequency weighted rank variability around rank frequency weighted mean
  #   sw_rel[i] <- sqrt(sum(pi*(ris-mw[i])^2)/sum(pi))
  #
  # }
  #out$ik
  #out$minmax
  #t <- table(out$ik)
  ####################################

  p<- max_pi[which(max_pi>=0.4)]
  names(p) <- which(max_pi>=0.4)
  p <- as.table(p)
  #p <- as.table(max_pi[ind])

  if(flag_show_fig){
    plot.new()
    frame()
    pdf("Figure_nri_check.pdf", width=6,height=5,paper='special')
    plot(p*100,ylim = c(0,100), xlab = "feature index",ylab = "frequency [%]",main = "Max. occupation frequency for RI/NRI features")
    #grid(nx = NA, ny = 4, col = "lightgray", lty = "dotted",
    #     lwd = par("lwd"), equilogs = TRUE)
    abline(h = 0, v = NA, col = "gray60")
    abline(h = 100, v = NA, col = "gray60")
    abline(h = 50, v = NA, col = "gray60",lty = "dashed")
    dev.off()
    print("Save figure to Figure_nri_check.pdf")
  }

  nri <- p*100

  #rownames(dat)[as.integer(names(p))]

  if(flag_show_fig){
    plot.new()
    frame()
    pdf("Figure_example_qn.pdf", width=6,height=5,paper='special')
    low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
    up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))
    if(length(ip)>1){
      boxplot(qn.dat,col=(c("gold")),notch=TRUE, xlab = "column", main = "Normalized data and maximum value")
      matlines(t(qn.dat[ip,]),type="l",col=c(4),ylim = c(low,up),xlab = "column",ylab = "normalized data")
      matlines(t(mbqn.dat[ip,]),type="l",col=c(3),ylim = c(low,up),xlab = "column",ylab = "normalized data")
    }else{
      plot(qn.dat[ip,],type="l",col=c(4),ylim = c(low,up),xlab = "column",ylab = "normalized data")
      boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = "column", main = "Normalized data and maximum value")
      lines(mbqn.dat[ip,],type="l",col=c(3),ylim = c(low,up),xlab = "column",ylab = "normalized data")
    }
    legend(x = "topright",legend=(c("qn-data","mbqn-data")),fill = c(4,3),bty= "n",cex=1,ncol=1)
    dev.off()
    print("Save figure to Figure_example_qn.pdf")
  }

  # mbqn.boxplot(qn.dat, irow = ip, vals = NULL, xlab = "column", ylab = "normalized data", main = "boxplot",filename = "Figure_example_qn2.df", type = "l")
  return(list(max_p = max_p,ip = ip,qn_dat = qn.dat, nri = nri, ind_var0 = ind_var0))

}
