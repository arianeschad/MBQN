#' Check data matrix for rank invariant (ri) /nearly rank invariant (nri) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param mean_fun median or mean, if left empty, quantile normalization
#' is applied without balancing the data
# #' @param qlow lower quantile
# #' @param qup upper quantile, default 1
#' @param flag_show_fig flag to specify whether results are plotted to figure, default = TRUE
#' @param low_thr lower threshold for NRI frequency, default = 0.2
#' @param filename string for naming figures, default = NULL
#' @param feature_index index of a feature of interest that is plotted in the boxplot; default NULL
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
#mbqn.check_saturation <- function(dat, FUN = mean_fun, qlow, qup, flag_show_fig = TRUE, low_thr = 0.2, filename = NULL){
mbqn.check_saturation <- function(dat, FUN = mean_fun, flag_show_fig = TRUE, low_thr = 0.2, filename = NULL,feature_index = NULL){

  # dat <- matrix(rnorm(20000,0,1),nrow = 2000, ncol = 10)
  # dat[1,] <- dat[1,]+4
  # dat[2,] <- dat[2,]+2
  # dat[3,] <- dat[3,]+1

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # quantile normalisation and its standard deviation
  qn.dat <- MBQN:::mbqn(x = dat,FUN = NULL)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  # quantile normalisation and its standard deviation with the limma function
  limma_qn.dat <- limma::normalizeBetweenArrays(object = dat ,method = "quantile")
  # limma::normalizeQuantiles()
  limma_s.qn <- apply(limma_qn.dat, 1, sd, na.rm=TRUE)


  # mean balanced quantile normalisation
  mbqn.dat <- MBQN:::mbqn(dat,FUN = FUN, na.rm = TRUE)

  ##############################################################################
  ## rank frequencies for each feature after QN (top-down)
  # assign NAs to 0 rank

  out <- MBQN:::get_kminmax(X = qn.dat,k = N, flag = "max")

  tdummy = lapply(
    1:N,
    function(i)
    {
      table(which(out$ik==i,arr.ind =T)[,1]*(!is.na(dat[i,])),useNA = 'ifany')
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

      # this part counts rank occupation frequencies if NAs are present
      #if(any(names(pi[[i]])==0)){
      #  Mp <- M - as.numeric(pi[[i]]["0"]*M)
      #  pi[[i]][ki] <- pi[[i]][ki]*M/Mp
      #}

      bla <- pi[[i]][which(pi[[i]]==max(pi[[i]][ki], na.rm =T))]
      bla
    }
  )

  max_pi_vals <- lapply(max_pi,function(i) unique(i))

  p <- max_pi_vals[which(max_pi_vals>=low_thr)]
  p <- unlist(p)
  names(p) <- as.character(which(max_pi_vals>=low_thr))
  p <- as.table(p)

  if(length(p)>1){
    warning('There might be multiple RI/NRI features!')
  }


  max_p <- max(p)*100
  ip <- as.integer(names(which(p*100==max_p)))

  freq_ismissing = sum(is.na(qn.dat[ip,]))/dim(qn.dat)[2] # cnt how often protein is missing

  print(paste('Maximum frequency of RI/NRI feature(s): ',max_p,"%"))

  #########################################################################################
  # which features have zero variation after QN

  ind_var0 <- which(s.qn==0)

  # how often are data present for these features
  not_nas <- apply(qn.dat,1,function(x) length(which(!is.na(x))))/M

  nri <- p*100

  p <- rbind(p,not_nas[as.numeric(names(p))])
  rownames(p) <- c("occupation.freq", "sample.coverage")

  not_nas <- not_nas[ind_var0]
  names(not_nas) <- ind_var0
  not_nas <- as.table(not_nas)


  ####### Graphical output #########

  # Occupation frequencies of RI/NRI features and sample coverage of zero variance features

  if(flag_show_fig){
    plot.new()
    frame()

    if(is.null(filename)) {
      fig1.name <- "Figure_nri_check.pdf"
    }else{fig1.name <- paste0("Figure_nri_check_", filename ,".pdf")}

    pdf(fig1.name, width=8,height=4,paper='special')

    plot(as.table(p[1,])*100, xlim = c(0,N), xaxt = "n", col = c(1),ylim = c(0,100), xlab = "feature index",
         ylab = "frequency [%]",cex.main=1,cex.lab=1,cex.lab =1, cex.axis=1,
         main = "Maxiumum occupation frequency of RI/NRI features\n & sample coverage of normalized features with zero variance ")

    lines(not_nas*100,xlim = c(0,N),ylim = c(0,100), lty=2, col = c(2))

    xticklabels <- c(names(not_nas),colnames(p))
    axis(1,at=as.integer(xticklabels), labels = FALSE)
    text(as.integer(xticklabels), par("usr")[3] - 5, labels = xticklabels, cex=0.75, srt = 90, pos = 1, xpd = TRUE)

    abline(h = 0, v = NA, col = "gray70")
    abline(h = 100, v = NA, col = "gray70")
    abline(h = 50, v = NA, col = "gray70")#,lty = "dashed")
    abline(h = low_thr*100, v = NA, col = c(4),lty = "dashed")
    op <- par(cex = 0.8)
    legend(x = "topright",legend=c("RI/NRI features","var=0 features","threshold"),col=c(1,2,4), lty=c(1:2,2), cex=0.8,bty = "n")

    dev.off()
    print(paste("Save figure to ",fig1.name))

    ###########################################################################################
    # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn

    plot.new()
    frame()
    if(is.null(filename)) {
      fig2.name <- "Figure_example_qn.pdf"
    }else{fig2.name <- paste0("Figure_example_qn_", filename ,".pdf")}

    pdf(fig2.name, width=8,height=4,paper='special')
    low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
    up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))
    if(length(ip)>1){
      boxplot(qn.dat,col=(c("gold")),notch=TRUE, xlab = "sample", main = "Quantile normalized data and maximum RI/NRI feature",cex.main = 0.8, outcex=0.3, xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      matlines(t(qn.dat[ip,]),type="l",col=c(4),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
      matlines(t(mbqn.dat[ip,]),type="l",col=c(2),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
    }else{
      plot(qn.dat[ip,],type="l",col=c(4),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = "sample", main = "Quantile normalized data and maximum RI/NRI feature",cex.main = 0.8, outcex=0.3,xaxt = "n")
      lines(mbqn.dat[ip,],type="l",col=c(2),ylim = c(low,up),xlab = "sample",ylab = "normalized data",xaxt = "n")
    }
    op <- par(cex = 0.8)
    legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),col=c(4,2), lty=1, cex=1,bty = "n") #fill = c(4,2),bty= "n",cex=1,ncol=1)
    dev.off()
    print(paste("Save figure to ",fig2.name))
    ###########################################################################################

    # boxplot of mbqn data and maximum RI/NRI feature after qn and mbqn

    plot.new()
    frame()
    if(is.null(filename)) {
      fig4.name <- "Figure_example_mbqn.pdf"
    }else{fig4.name <- paste0("Figure_example_mbqn_", filename ,".pdf")}

    pdf(fig4.name, width=8,height=4,paper='special')
    low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
    up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))
    if(length(ip)>1){
      boxplot(mbqn.dat,col=(c("gold")),notch=TRUE, xlab = "sample", main = "MBQN data and maximum RI/NRI feature",cex.main = 0.8, outcex=0.3, xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      matlines(t(qn.dat[ip,]),type="l",col=c(4),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
      matlines(t(mbqn.dat[ip,]),type="l",col=c(2),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
    }else{
      plot(qn.dat[ip,],type="l",col=c(4),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      boxplot(mbqn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = "sample", main = "Quantile normalized data and maximum RI/NRI feature",cex.main = 0.8, outcex=0.3,xaxt = "n")
      lines(mbqn.dat[ip,],type="l",col=c(2),ylim = c(low,up),xlab = "sample",ylab = "normalized data",xaxt = "n")
    }
    op <- par(cex = 0.8)
    legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),col=c(4,2), lty=1, cex=1,bty = "n") #fill = c(4,2),bty= "n",cex=1,ncol=1)
    dev.off()
    print(paste("Save figure to ",fig4.name))

    ###########################################################################################

    plot.new()
    frame()

    if(is.null(filename)) {
      fig3.name <- "Figure_qn_nri_check.pdf"
    }else{fig3.name <- paste0("Figure_qn_nri_check_", filename ,".pdf")}

    pdf(fig3.name, width=10,height=9,paper='special')
    par(mfrow=c(2,1))

    # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn

    low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
    up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))
    if(!is.null(feature_index)) ip <- unique(c(ip,feature_index))

    if(length(ip)>1){
      boxplot(qn.dat,col=(c("gold")),notch=TRUE, ylim = c(low,up),xlab = "sample", main = "Quantile normalized data and maximum RI/NRI feature",cex.main = 0.8, outcex=0.4, xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .9)
      matlines(t(qn.dat[ip,]),type="l",col=c(4),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
      matlines(t(mbqn.dat[ip,]),type="l",col=c(2),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
    }else{
      plot(qn.dat[ip,],type="l",col=c(4),ylim = c(low,up),xlab = "sample",ylab = "normalized data", xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .9)
      boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = "sample", main = "Quantile normalized data and maximum RI/NRI feature",cex.main = 0.8, outcex=0.4,xaxt = "n")
      lines(mbqn.dat[ip,],type="l",col=c(2),ylim = c(low,up),xlab = "sample",ylab = "normalized data",xaxt = "n")
    }
    par(cex.lab=.8) # y-axis
    par(cex.axis=.8) # x-axis
    #op <- par(cex = 0.8)
    legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),col=c(4,2), lty=1, cex=0.8, bty = "n")
    fig_label("A.", cex=1.7)

    # ri/nri check and zero variance qn features

    plot(as.table(p[1,])*100, xlim = c(0,N), xaxt = "n", col = c(1),ylim = c(0,100), xlab = "feature index",
         ylab = "frequency [%]",cex.main = .85, cex.lab = 1,cex.lab = 1, cex.axis = 1,
         main = "Maxiumum occupation frequency of RI/NRI features\n & sample coverage of normalized features with zero variance ")

    lines(not_nas*100,xlim = c(0,N),ylim = c(0,100), lty=2, col = c(2))

    xticklabels <- c(names(not_nas),colnames(p))
    axis(1,at = as.integer(xticklabels), labels = FALSE, cex.axis = .8)
    text(as.integer(xticklabels), par("usr")[3] - 5, labels = xticklabels, cex=0.75, srt = 90, pos = 1, xpd = TRUE)

    abline(h = 0, v = NA, col = "gray70")
    abline(h = 100, v = NA, col = "gray70")
    abline(h = 50, v = NA, col = "gray70")#,lty = "dashed")
    abline(h = low_thr*100, v = NA, col = c(4),lty = "dashed")
    op <- par(cex = .8)
    legend(x = "topright",legend=c("RI/NRI features","var=0 features","threshold"),col=c(1,2,4), lty=c(1:2,2), cex=0.8, bty = "n")

    fig_label("B.", cex=2)

    dev.off()
    print(paste("Save figure to ",fig3.name))

  }

  return(list(max_p = max_p, ip = ip, nri = nri, var0_feature = ind_var0))

}

# helper function from
# https://logfc.wordpress.com/2017/03/15/adding-figure-labels-a-b-c-in-the-top-left-corner-of-the-plotting-region/

fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {

  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright",
                          "left", "center", "right",
                          "bottomleft", "bottom", "bottomright"))

  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")

    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    }
  }

  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }

  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100

  x1 <- switch(pos,
               topleft     =x[1] + sw,
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)

  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)

  old.par <- par(xpd=NA)
  on.exit(par(old.par))

  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}

