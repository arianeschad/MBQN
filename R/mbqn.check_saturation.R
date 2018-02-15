#' Test data for saturated upper tail
# quantile-normalization of data benefits from mean-balancing
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param mean_fun mean or median, if left empty, quantile normalization
#' is applied without balancing the data
#' @param qlow lower quantile
#' @param qup upper quantile, default 1
#' @param flag_show_fig flag to specify whether plot results to figure, default = TRUE
#' @inheritParams mbqn
#' @details Rank data and test if lower and upper intensity tails are
#' dominated by few feature types/proteins. Compute a quantile
#' normalization without and with mean-balancing and check standard
#' deviation of normalized data entries located in the tails
#' @return A matrix of mean- or median-balanced quantile normalized data
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile
#' normalization for processing label-free quantitative proteomics
#' data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples mbqn.check_saturation(dat, mean)
#' @description Test if few rows of the data matrix dominate the upper tail. This script uses normalize.quantiles() from the
#' package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
#' @author A. Schad, \email{ariane.schad@zbsa.de}


# dat <- matrix(rnorm(20000,0,1),nrow = 2000, ncol = 10)
# dat[1,] <- dat[1,]+4
# dat[2,] <- dat[2,]+2
# dat[3,] <- dat[3,]+1

mbqn.check_saturation <- function(dat, FUN = mean_fun, qlow, qup, flag_show_fig = TRUE){

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # row means
  m.dat <- apply(dat,1,FUN, na.rm=TRUE)

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  #mean balanced quantile normalisation and its standard deviation
  mbqn.dat <- mbqn(dat,FUN = FUN, na.rm = TRUE)
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
  s<- matrix(NA,nrow = 20,nco=dim(ir)[1])
  for (i in 1:20){
    test <- sapply(1:dim(ir)[2],function(i) sample(ir[,i],dim(ir)[1],replace = FALSE))
    s[i,] <- apply(test, 1, sd, na.rm=TRUE)
  }
  s.m <- colMeans(s)
  s.s <- apply(s,2,sd,na.rm = TRUE)

  plot(s.ir,type="l",col=c(2))
  lines(s.m)

  ## find top 5 of largest values in each column
  out <- get_kminmax(X = dat,k = 5, flag = "max")
  out$ik
  out$minmax
  t <- table(out$ik)
  ####

  # quantile normalisation
  qn.dat <- mbqn(x = dat,FUN = NULL)

  mbqn.dat <- mbqn(dat,FUN = FUN, na.rm = TRUE)
  out=get_kminmax(X = qn.dat,k = 5,flag = "max")
  t = table(out$ik)
  p <- t/M # frequency of the most frequent highest-expressed protein
  max_p <- max(p)
  ip <- which(t==max(t)) # most frequent top-level protein

  if(length(ip)>1){
    #ip = ip[1]
    warning('There are multiple high level proteins!')
    flag_ambigious<-TRUE
  }
  freq_ismissing = sum(is.na(qn.dat[ip,]))/dim(qn.dat)[2] # cnt how often protein is missing

  print(paste('Frequency of single highest protein: ',max_p*100," %"))

 if(flag_show_fig){
  pdf("Figure_example_qn.pdf", width=6,height=5,paper='special')
  plot(qn.dat[ip,],type="l",col=c(4),ylim = c(-4.5,6),xlab = "column",ylab = "QN data")
  boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = "column", main = "Normalized data and maximum value")
  lines(mbqn.dat[ip,],type="l",col=c(3),ylim = c(-4.5,6),xlab = "column",ylab = "QN data")
  legend(x = "topright",legend=(c("qn-data","mbqn-data")),fill = c(4,3),bty= "n",cex=1,ncol=1)
  dev.off()
  print("Save figure to Figure_example_qn.pdf")
}

  return(list(max_p =max_p,ip = ip,qn_dat = qn_dat))
}
