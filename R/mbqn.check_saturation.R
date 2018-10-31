#' Check data matrix for rank invariant (ri) /nearly rank invariant (nri) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param FUN median or mean, if left empty, quantile normalization
#' is applied without balancing the data
# #' @param qlow lower quantile
# #' @param qup upper quantile, default 1
#' @param flag_show_fig flag to specify whether results are plotted to figure, default = TRUE
#' @param low_thr lower threshold for NRI frequency, default = 0.2
#' @param filename string for naming figures, default = NULL
#' @param feature_index index of a feature of interest that is plotted in the boxplot; default NULL
#' @param method Package used to compute quantile normalization; default NULL - use the preprocessCore package ; "limma" uses the Limma package
#' @inheritParams mbqn
#' @details Rank data and check if lower and upper intensity tails are
#' dominated by few feature. Compute a quantile
#' normalization without and with mean-balancing and check standard
#' deviation of normalized data entries located in the tails
#' @return A matrix of median- or mean-balanced quantile normalized data
#' @keywords quantile normalization, proteomics
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
#' @examples mbqn.check_saturation(dat, mean)
#' @description Test if few rows of the data matrix dominate the upper tail. This script uses normalize.quantiles() from the
#' package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' 2017
#' @export
mbqn.check_saturation <- function(dat, FUN = NULL,
                                  flag_show_fig = TRUE,
                                  low_thr = 0.2,
                                  filename = NULL,
                                  feature_index = NULL,
                                  method = NULL){

  # if FUN is not specified, use median!
  if(is.null(FUN)) FUN <- median

  N <- dim(dat)[1] #number of rows
  M <- dim(dat)[2] #number of cols

  # plot options
  cex.main <- 1.2
  cex.lab <- 1.
  cex.legend <- 1.
  cex.axis <- 1.

  # quantile normalisation and its standard deviation
  qn.dat <- MBQN:::mbqn(x = dat,FUN = NULL)
  s.qn <- apply(qn.dat, 1, sd, na.rm=TRUE)

  # quantile normalisation and its standard deviation computed with the limma function
  # The results from limma::normalizeBetweenArrays normalization are equal, up to numerical differences,
  # to that of preprocessCore::normalize.quantiles!

  # mean balanced quantile normalisation, optionally with limma qn function
  mbqn.dat <- MBQN:::mbqn(dat,FUN = FUN, na.rm = TRUE, method = method)

  ##############################################################################
  ## Rank frequencies for each feature after QN (top-down)
  # & assign NAs to 0 rank

  out <- MBQN::get_kminmax(X = qn.dat,k = N, flag = "max")

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
      bla <- pi[[i]][which(pi[[i]]==max(pi[[i]][ki], na.rm =T))]
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

    #    xlimit <- c(max(as.numeric(names(which.min(nri)))-3,0),as.numeric(names(which.max(nri)))+3)
    xlimit <- c(min(as.numeric(names(nri)))-3,max(as.numeric(names(nri)))+3)

    plot(as.table(p[1,])*100,
         #xlim = c(0,N),
         xlim = xlimit,
         xaxt = "n",
         col = c(1),ylim = c(0,100),
         xlab = "feature index",
         ylab = "frequency [%]",
         cex.main=cex.main,cex.lab =cex.lab ,cex.axis=cex.axis,
         main = "Maxiumum occupation frequency of RI/NRI features\n & sample coverage of normalized features with zero variance ")

    lines(not_nas*100,xlim = c(0,N),ylim = c(0,100), lty=2, col = c(2))

    xticklabels <- c(names(not_nas),colnames(p))
    axis(1,at=as.integer(xticklabels), labels = FALSE)
    text(as.integer(xticklabels), par("usr")[3] - 5, labels = xticklabels,
         cex=0.75, srt = 90, pos = 1, xpd = TRUE)

    abline(h = 0, v = NA, col = "gray70")
    abline(h = 100, v = NA, col = "gray70")
    abline(h = 50, v = NA, col = "gray70")#,lty = "dashed")
    abline(h = low_thr*100, v = NA, col = c(4),lty = "dashed")
    op <- par(cex = 0.8)
    legend(x = "topleft",legend=c("RI/NRI features","var=0 features","threshold"),
           col=c(1,2,4), lty=c(1:2,2), cex=cex.legend ,bty = "n")

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
      boxplot(qn.dat,col=(c("gold")),notch=F, xlab = "sample",
              main = "Quantile normalized data and maximum RI/NRI feature",
              cex.main = cex.main, outcex=0.3, xaxt = "n")

      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      matlines(t(qn.dat[ip,]),type="b",col=c(4),ylim = c(low,up),
               xlab = "sample",ylab = "normalized data", xaxt = "n")
      matlines(t(mbqn.dat[ip,]),type="b",col=c(2),ylim = c(low,up),
               xlab = "sample",ylab = "normalized data", xaxt = "n")
    }else{
      plot(qn.dat[ip,],type="b",col=c(4),ylim = c(low,up),xlab = "sample",
           ylab = "normalized data", xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=F, xlab = "sample",
              main = "Quantile normalized data and maximum RI/NRI feature",
              cex.main = cex.main, outcex=0.3,xaxt = "n")
      lines(mbqn.dat[ip,],type="b",col=c(2),ylim = c(low,up),xlab = "sample",
            ylab = "normalized data",xaxt = "n")
    }
    op <- par(cex = 0.8)
    legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),
           col=c(4,2), lty=1, cex=cex.legend ,bty = "n") #fill = c(4,2),bty= "n",cex=1,ncol=1)
    dev.off()
    print(paste("Save figure to ",fig2.name))
    ###########################################################################################

    # boxplot of mbqn data and maximum RI/NRI feature after qn and mbqn

    plot.new()
    frame()
    par(mfrow = c(1,1), cex.lab = 1.5)

    if(is.null(filename)) {
      fig4.name <- "Figure_example_mbqn.pdf"
    }else{fig4.name <- paste0("Figure_example_mbqn_", filename ,".pdf")}

    pdf(fig4.name, width=8,height=4,paper='special')

    low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
    up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))

    if(length(ip)>1){
      boxplot(mbqn.dat,col=(c("gold")),notch=F, xlab = "sample",
              ylab = "normalized intensity",
              main = "MBQN data and maximum RI/NRI feature",
              cex.main = cex.main,outcex=0.2, xaxt = "n")
      #par(cex.lab = 1.5)
      #par(cex.main = 0.5)
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = 1.2)

      matlines(t(qn.dat[ip,]),type="b",col=c(4),ylim = c(low,up),
               xlab = "sample",ylab = "normalized intensity", xaxt = "n")
      matlines(t(mbqn.dat[ip,]),type="b",col=c(2),ylim = c(low,up),
               xlab = "sample",ylab = "normalized intensity", xaxt = "n")
      par(cex.lab=cex.lab) # is for y-axis
      par(cex.axis=cex.lab) # is for x-axis

    }else{
      plot(qn.dat[ip,],type="b",col=c(4),ylim = c(low,up),
           #ylim = c(35,36),
           xlab = "sample",
           ylab = "normalized data", xaxt = "n",
           cex.lab = cex.lab,
           lwd  =1.5)
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = 1)
      boxplot(mbqn.dat,col=(c("gold")),add = TRUE,notch=F,
              # xlab = "sample",
              main = "Median-balanced quantile normalized data\n & maximum RI/NRI features",
              cex.main = cex.main, outcex=0.3,xaxt = "n")
      #lines(mbqn.dat[ip,],type="l",col=c(2),ylim = c(low,up)) # ,
      #     xlab = "sample",ylab = "normalized data",xaxt = "n")

      lines(mbqn.dat[ip,],col=c(2),ylim = c(low,up), type="b", lwd=1.5)#, #type="l"
      # lty=1, pch=1)

    }
    op <- par(cex = 0.8)
    legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),
           col=c(4,2), lty=1, lwd = 1.5, cex=cex.legend, bty = "n") #fill = c(4,2),bty= "n",cex=1,ncol=1)
    dev.off()
    print(paste("Save figure to ",fig4.name))

    ###########################################################################################

    plot.new()
    frame()

    if(is.null(filename)) {
      fig3.name <- "Figure_qn_nri_check.pdf"
    }else{fig3.name <- paste0("Figure_qn_nri_check_", filename ,".pdf")}

    pdf(fig3.name, width=8,height=7,paper='special')
    par(mfrow=c(2,1))

    # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn

    low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
    up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))
    if(!is.null(feature_index)) ip <- unique(c(ip,feature_index))

    if(length(ip)>1){
      boxplot(qn.dat,col=(c("gold")),notch=F, ylim = c(low,up),
              xlab = "sample",
              main = "Quantile normalized data and maximum RI/NRI feature",
              cex.main = cex.main, outcex=0.4, xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .9)
      matlines(t(qn.dat[ip,]),type="b",pch=1, col=c(4),ylim = c(low,up),
               xlab = "sample",ylab = "normalized data", xaxt = "n")
      matlines(t(mbqn.dat[ip,]),type="b",pch=1,col=c(2),ylim = c(low,up),
               xlab = "sample",ylab = "normalized data", xaxt = "n")
    }else{
      plot(qn.dat[ip,],col=c(4), type="b", ylim = c(low,up),xlab = "sample",
           ylab = "normalized data", xaxt = "n")
      axis(1, at = c(1:M), labels = c(1:M), cex.axis = .9)
      boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=F,
              xlab = "sample", main = "Quantile normalized data and maximum RI/NRI feature",
              cex.main = cex.main, outcex=0.4,xaxt = "n")
      #lines(mbqn.dat[ip,],col=c(2),ylim = c(low,up), type="l", xaxt = "n")
      lines(mbqn.dat[ip,],col=c(2),ylim = c(low,up), type="b", lwd=1.5, xaxt = "n")#, #type="l"
      # lty=1, pch=1)

    }
    par(cex.lab=.8) # y-axis
    par(cex.axis=.8) # x-axis
    #op <- par(cex = 0.8)
    legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),col=c(4,2),
           lty=1, cex=0.8, bty = "n")
    fig_label("A.", cex=1.7)

    # ri/nri check and zero variance qn features
    plot(as.table(p[1,])*100, xlim = xlimit, #xlim = c(0,N),
         xaxt = "n", col = c(1),
         ylim = c(0,100), xlab = "feature index",
         ylab = "frequency [%]",cex.main = .85,
         cex.lab = cex.lab, cex.axis = 1,cex.main =cex.main,
         main = "Maximum occupation frequency of RI/NRI features\n & sample coverage of normalized features with zero variance ")

    lines(not_nas*100,xlim = c(0,N),ylim = c(0,100), lty=2, col = c(2))

    xticklabels <- c(names(not_nas),colnames(p))
    axis(1,at = as.integer(xticklabels), labels = FALSE, cex.axis = .8)
    text(as.integer(xticklabels), par("usr")[3] - 5, labels = xticklabels,
         cex=0.75, srt = 90, pos = 1, xpd = TRUE)

    abline(h = 0, v = NA, col = "gray70")
    abline(h = 100, v = NA, col = "gray70")
    abline(h = 50, v = NA, col = "gray70")#,lty = "dashed")
    abline(h = low_thr*100, v = NA, col = c(4),lty = "dashed")
    op <- par(cex = .8)
    legend(x = "topleft",legend=c("RI/NRI features","var=0 features","threshold"),
           col=c(1,2,4), lty=c(1:2,2), cex=0.8, bty = "n")

    fig_label("B.", cex=2)

    dev.off()
    print(paste("Save figure to ",fig3.name))

  }

  return(list(max_p = max_p, ip = ip, nri = nri, var0_feature = ind_var0))

}

#####################################################################
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

