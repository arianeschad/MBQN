#' Check data matrix for rank invariant (ri) /nearly rank invariant (nri) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param FUN median or mean, if left empty, quantile normalization
#' is applied without balancing the data
# #' @param qlow lower quantile
# #' @param qup upper quantile, default 1
#' @param flag_show_fig flag to specify whether results are plotted to figure, default = TRUE
#' @inheritParams mbqn
#' @details Rank data and check if lower and upper intensity tails are
#' dominated by few feature. Compute a quantile
#' normalization without and with mean-balancing and check standard
#' deviation of normalized data entries located in the tails
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
mbqn.check_saturation_v0 <- function(dat, FUN = NULL, flag_show_fig = TRUE){
#qlow = NULL, qup = NULL,

  # if FUN is not specified, use median!
  if(is.null(FUN)) FUN <- median

  # dat <- matrix(rnorm(20000,0,1),nrow = 2000, ncol = 10)
  # dat[1,] <- dat[1,]+4
  # dat[2,] <- dat[2,]+2
  # dat[3,] <- dat[3,]+1

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # row means
  m.dat <- apply(dat,1,FUN, na.rm=TRUE)

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  #mean balanced quantile normalisation and its standard deviation
  mbqn.dat <- MBQN:::mbqn(dat,FUN = FUN, na.rm = TRUE)
  s.mbqn <- apply(mbqn.dat, 1, sd, na.rm=TRUE)

  # compute (tied) ranks
  ir <- sapply(1:M,function(i) rank(dat[,i], na.last = TRUE, ties.method = c("average")))
  s.ir <- apply(ir, 1, sd, na.rm=TRUE)
  which(s.ir == 0)
  plot(s.ir,type = "l")
  p <- 5

  # find indices of the top n
  top <- sapply(1:M,function(i) which(ir[,i] %in% c((max(ir)-p):max(ir))))

  # Trick: shuffle ranks of each column and compute statistics
  s <- matrix(NA,nrow = 20,nco=dim(ir)[1])
  for (i in 1:20){
    test <- sapply(1:dim(ir)[2],function(i) sample(ir[,i],dim(ir)[1],replace = FALSE))
    s[i,] <- apply(test, 1, sd, na.rm=TRUE)
  }
  s.m <- colMeans(s)
  s.s <- apply(s,2,sd,na.rm = TRUE)

  plot(s.ir,type="l",col=c(2))
  lines(s.m)

  ## find top 5 of largest values in each column
  out <- MBQN:::get_kminmax(X = dat,k = 5, flag = "max")
  out$ik
  out$minmax
  t <- table(out$ik)
  ####

  # quantile normalisation
  qn.dat <- MBQN:::mbqn(x = dat,FUN = NULL)
  # mean balanced quantile normalisation
  mbqn.dat <- MBQN:::mbqn(dat,FUN = FUN, na.rm = TRUE)

  out=MBQN:::get_kminmax(X = qn.dat,k = 10,flag = "max")
  t = table(out$ik)
  p <- t/M # frequency of the most frequent highest-expressed protein
  max_p <- max(p)
  ip <- as.integer(names(which(t==max(t)))) # index of most frequent top-level protein

  if(length(ip)>1){
    #ip = ip[1]
    warning('There are multiple RI/NRI features!')
    flag_ambigious<-TRUE
  }
  freq_ismissing = sum(is.na(qn.dat[ip,]))/dim(qn.dat)[2] # cnt how often protein is missing

  print(paste('Maximum frequency of RI/NRI feature(s): ',max_p*100,"%"))

  if(flag_show_fig){
    plot.new()
    frame()
    pdf("Figure_nri_check.pdf", width=6,height=5,paper='special')
    plot(p*100,ylim = c(0,100), xlab = "feature index",ylab = "frequency [%]",main = "Occupation frequency for RI/NRI features")
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
      boxplot(qn.dat,col=(c("gold")),notch=F, xlab = "column", main = "Normalized data and maximum value")
      matlines(t(qn.dat[ip,]),type="l",col=c(4),ylim = c(low,up),xlab = "column",ylab = "normalized data")
      matlines(t(mbqn.dat[ip,]),type="l",col=c(3),ylim = c(low,up),xlab = "column",ylab = "normalized data")
    }else{
      plot(qn.dat[ip,],type="l",col=c(4),ylim = c(low,up),xlab = "column",ylab = "normalized data")
      boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=F, xlab = "column", main = "Normalized data and maximum value")
      lines(mbqn.dat[ip,],type="l",col=c(3),ylim = c(low,up),xlab = "column",ylab = "normalized data")
    }
    legend(x = "topright",legend=(c("qn-data","mbqn-data")),fill = c(4,3),bty= "n",cex=1,ncol=1)
    dev.off()
    print("Save figure to Figure_example_qn.pdf")
  }

  ##################################################################################
  # which features show no variation after QN
  ind_var0 <- which(apply(qn.dat,1,sd, na.rm=T)==0)
  # which of them contain mainly na's


  #mbqn.boxplot(qn.dat, irow = ip, vals = NULL, xlab = "column", ylab = "normalized data", main = "boxplot",filename = "Figure_example_qn2.df", type = "l")

  return(list(max_p =max_p,ip = ip,qn_dat = qn.dat, nri = nri, ind_var0 = ind_var0))
}
