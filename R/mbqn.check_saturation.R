#' Check data matrix for rank invariant (RI) and nearly rank invariant (NRI) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param FUN median, mean, or another function used to balance features across colums. If left empty, quantile normalization
#' is applied without balancing the data
# #' @param qlow lower quantile
# #' @param qup upper quantile, default 1
#' @param low_thr Numerical value for the lower threshold for NRI frequency, default = 0.5
#' @param feature_index Integer that indicates the index of a feature of interest that is plotted in the boxplot; default NULL
#' @param method Packagename containing function used to compute quantile normalization; default NULL - use the preprocessCore package ; "limma" uses the Limma package
#' @param show_fig Logical flag indicating whether results should be displayed; default = TRUE.
#' @param save_fig Logical to save figures to file
#' @param show_nri_only Logical to print and save only the RI/NRI detection graph; default FALSE
#' @param filename String for naming figures, default = NULL
#' @param verbose Logical for running function quiet
#' @param ... Optional plot arguments passed to \code{mbqn.boxplot}
#' @inheritParams mbqn
#' @inheritParams mbqn.boxplot
#' @importFrom grDevices dev.copy2pdf dev.off dev.size pdf
#' @details Rank data and check if lower and upper intensity tails are
#' dominated by few features. Apply quantile
#' normalization without and with mean-balancing and check the standard
#' deviation of normalized features located in the tails.
#' @return List with the components:
#' \item{\code{p}}{a matrix of the rank invariance frequencies and the sample coverage over all RI/NRI features}
#' \item{\code{max_p}}{value of the maximum rank invariance frequency in percent}
#' \item{\code{ip}}{index of the feature with maximum rank invariance frequency}
#' \item{\code{nri}}{a table of the rank invariance frequencies in percent of NRI/RI features}
#' \item{\code{var0_feature}}{index of features with zero sample variance after QN.}
# detected RI and NRI feature indices together with their rank invariance
# frequency and sample coverage, indices of features with zero variation after QN, and the feature index of maximum rank invariance frequency. The
# sample coverage measures 1-(relative amount of missing values) of the features. Optionally, a graphic
# depicts the detected RI/NRI features, frequencies, and sample coverage and boxplots of normalized data with RI/NRI features.
# #' @keywords quantile normalization, proteomics
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
#' @examples ## Check data matrix for RI and NRI features
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqn.check_saturation(X, mean, low_thr = 0.5, save_fig = FALSE)
#' @description Check data matrix for rows with intensity features which dominate the upper tail, i.e. for
#' features that have constant or nearly constant rank across samples.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
# 2017
#' @export mbqn.check_saturation
mbqn.check_saturation <- function(dat, FUN = NULL,
                                  low_thr = 0.5,
                                  feature_index = NULL,
                                  method = NULL,
                                  show_fig = TRUE,
                                  save_fig = TRUE,
                                  show_nri_only = FALSE,
                                  filename = NULL,verbose = TRUE,...){


  res  <- mbqn.get_nri_features(dat, FUN = FUN,
                                low_thr = low_thr,
                                method = method,
                                verbose = verbose)

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL, method = method, verbose = FALSE)
  mbqn.dat <- mbqn.nri(x = dat, FUN = median, low_thr = low_thr, method = method, verbose = FALSE)

  ####### Graphical output #########

  # plot options
  cex.main <- 1.2
  cex.lab <- 1.
  cex.legend <- .8
  cex.axis <- 1.

  # Occupation or rank invariance frequencies
  # and sample coverage of RI/NRI features
  if(show_fig){
    if(!is.null(res$nri)){

      current.dir = getwd()

      plot.new()
      frame()
      par(mfrow=c(1,1), mar = c(4,4,3,2))

      if(is.null(filename)) {
        fig1.name <- "Figure_nri_check.pdf"
      }else{fig1.name <- paste0("Figure_nri_check_", filename ,".pdf")}

      ylim <- c(min(as.numeric(names(res$nri)))-1,max(as.numeric(names(res$nri)))+2)

      dummy <- data.frame(frequency <- as.numeric(res$p[1,]),
                          feature <- as.integer(colnames(res$p)))
      colnames(dummy) <- c("frequency","feature index")

      layout(matrix(c(1,2), nrow = 1,ncol = 2),widths = c(3,2))
      par(mar = c(6,4,4,0), las = 2)
      plot(NA, type="n", xlim=c(0.5,100),
           ylim=ylim,
           bty='L',
           xlab="frequency [%]",
           ylab="feature index",
           yaxt="n", xaxt="n", xaxs="i")
      axis(side = 1, seq(0,100,10), las =1,cex.axis = 0.8)
      rect(xleft=0, ybottom=dummy$`feature index`-0.05,
           xright=dummy$frequency*100,
           ytop=dummy$`feature index`+0.05, col = 1)
      legend.txt <- "RI/NRI feature"
      abline(v = 100, h = NA, col = "black")
      abline(v = low_thr*100, h = NA, col = c(4),lty = "dashed")
      axis(side = 2, at = dummy$`feature index`,cex.axis = 0.8)
      legend.txt <- cbind(legend.txt, "threshold")
      ind <- NULL
      # zero variance features
      #if(length(res$var0_feature)>0)
      #  ind <- which(dummy$`feature index`==res$var0_feature)
      ind <- as.numeric(names(res$nri))
      if(length(ind)>0) {
        if(length(ind)==1) {
        coverage <- sum(!is.na(qn.dat[ind,]))/dim(qn.dat)[2]
        } else { coverage <- apply(!is.na(qn.dat[ind,]), 1,sum)/dim(qn.dat)[2]
        }
      }

      # coverage <- sum(as.numeric(!is.na(qn.dat[dummy$`feature index`[ind],])))/dim(qn.dat)[2]
      #rect(xleft=0, ybottom=dummy$`feature index`-0.05,
      #     xright=coverage*100, ytop=dummy$`feature index`+0.05,
      #     border = 2, col =2)
      #  text(rep(95,length(ind)), dummy$`feature index`+0.5,
      #       labels = paste0(as.character(coverage*100),"%"),
      #       col="red" ,cex = 0.8)
      #rect(xleft=0, ybottom=dummy$`feature index`[ind]-0.05,
      #     xright=dummy$frequency[ind]*100, ytop=dummy$`feature index`[ind]+0.05,
      #     border = 2, col =2)

      par(mar = c(6,0.1,4,0))
      plot(NA, type="n", xlim=c(0,100),
           ylim=ylim,
           bty='n',
           xlab="",
           ylab="",
           yaxt="n",
           xaxt="n")
      if(length(ind)>0){
        text(rep(5,length(ind)), dummy$`feature index`,
             labels = paste0(as.character(coverage*100),"%"),
             col="red" ,cex = 0.8)
        legend.txt <- cbind(legend.txt, "feature coverage")
      }

      legend(x = 10, y = ylim[2],
             #x = "topleft",
             xpd = T,
             #legend=c("RI/NRI feature","feature coverage","threshold"),
             legend = legend.txt,
             col=c(1,4,2), lty=c(1,2,1), cex=cex.legend ,bty = "n", y.intersp = 0.3)
      if(save_fig){
        dev.copy2pdf(file=file.path(current.dir,fig1.name),width=6,height=11,out.type = "pdf", paper="a4")
        if(verbose) print(paste("Save figure to ",fig1.name))
      }

      ############################################
      #   plot(as.table(res$p[1,])*100,
      #      #xlim = c(0,N),
      #      xlim = xlimit,
      #      xaxt = "n",
      #      col = c(1),ylim = c(0,100),
      #      xlab = "feature index",
      #      ylab = "frequency [%]",
      #      cex.main=cex.main,cex.lab =cex.lab ,cex.axis=cex.axis,
      #      main = "Maxiumum occupation frequency of RI/NRI features\n & sample coverage of normalized features with zero variance ")
      #
      # lines(not_nas*100,xlim = c(0,N),ylim = c(0,100), lty=1, col = c(2))
      #
      # xticklabels <- c(names(not_nas),colnames(p))
      # axis(1,at=as.integer(xticklabels), labels = FALSE)
      # text(as.integer(xticklabels), par("usr")[3] - 5, labels = xticklabels,
      #      cex=0.75, srt = 90, pos = 1, xpd = TRUE)
      #
      # abline(h = 0, v = NA, col = "gray70")
      # abline(h = 100, v = NA, col = "gray70")
      # abline(h = 50, v = NA, col = "gray70")#,lty = "dashed")
      # abline(h = low_thr*100, v = NA, col = c(4),lty = "dashed")
      # op <- par(cex = 0.8)
      # legend(x = "topleft",inset=c(0,0.01),xpd = T, legend=c("RI/NRI features","coverage of RI features","threshold"),
      #        col=c(1,2,4), lty=c(1,1,2), cex=cex.legend ,bty = "n")
      #
      # if(save_fig){
      #   dev.copy2pdf(file=file.path(current.dir,fig1.name),width=8,height=4,out.type = "pdf", paper="a4r")
      #   if(verbose) print(paste("Save figure to ",fig1.name))
      # }
      ###########################################################################################
      # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn
      if(!show_nri_only){
      dev.off()
      plot.new()
      frame()
      #par(mfrow=c(1,1), mar = c(4,4,3,2))


      low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
      up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))

      mbqn.boxplot(qn.dat,
                   vals = data.frame(QN.feature = qn.dat[res$ip,], MBQN.feature = mbqn.dat[res$ip,]),
                   ylim = c(low,up),
                   ylab = "normalized intensity",
                   main = "QN data with unbalanced \n and balanced maximum RI/NRI feature",
                   cex.main = 0.8,...)

      # if(length(ip)>1){
      #   boxplot(qn.dat,col=(c("gold")),notch=F, xlab = "sample",
      #           main = "Quantile normalized data and maximum RI/NRI feature",
      #           cex.main = cex.main, outcex=0.3, ylab = "normalized data", xaxt = "n")
      #
      #   axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      #   matlines(t(qn.dat[ip,]),type="b",col=c(4),ylim = c(low,up), xaxt = "n")
      #   matlines(t(mbqn.dat[ip,]),type="b",col=c(2),ylim = c(low,up), xaxt = "n")
      # }else{
      #     #lines(mbqn.dat[ip,],type="b",col=c(2),ylim = c(low,up),xaxt = "n")

      # plot(qn.dat[ip,],type="b",col=c(4),ylim = c(low,up),xlab = "sample",
      #       ylab = "normalized data", xaxt = "n", cex.lab = 1.2,
      #       main = "Quantile normalized data and maximum RI/NRI feature",
      #       cex.main = cex.main)
      #  axis(1, at = c(1:M), labels = c(1:M), cex.axis = .8)
      #  boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=F, outcex=0.3,xaxt = "n")
      #  lines(mbqn.dat[ip,],type="b",col=c(2),ylim = c(low,up),xaxt = "n")
      #}
      # op <- par(cex = 0.8)
      # legend(x = "bottomright",legend=(c("QN feature","MBQN feature")),
      #         col=c(4,2), lty=1, cex=cex.legend ,bty = "n") #fill = c(4,2),bty= "n",cex=1,ncol=1)

      if(is.null(filename)) {
        fig2.name <- "Figure_example_qn.pdf"
      }else{fig2.name <- paste0("Figure_example_qn_", filename ,".pdf")}

      if(save_fig){
        dev.copy2pdf(file=file.path(current.dir,fig2.name),width=10,height=5,out.type = "pdf", paper="a4r")
        if(verbose) print(paste("Save figure to ",fig2.name))
      }
      ###########################################################################################

      # boxplot of qn data and with balanced qn RI/NRI features
      plot.new()
      frame()

      low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
      up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))

      mbqn.boxplot(mbqn.dat,
                   vals = data.frame(QN.feature = qn.dat[res$ip,], MBQN.feature = mbqn.dat[res$ip,]),
                   ylim = c(low,up),
                   ylab = "normalized intensity",
                   main = "MBQN data with unbalanced \n and balanced maximum RI/NRI feature",
                   cex.main = 0.8,...)

      if(is.null(filename)) {
        fig4.name <- "Figure_example_mbqn.pdf"
      }else{fig4.name <- paste0("Figure_example_mbqn_", filename ,".pdf")}

      # par(mfrow = c(1,1), cex.lab = 1.5)
      # if(length(ip)>1){
      #   boxplot(mbqn.dat,col=(c("gold")),notch=FALSE, xlab = "sample",
      #           ylab = "normalized intensity",
      #           main = "MBQN data and maximum RI/NRI feature",
      #           cex.main = cex.main,outcex=0.2, xaxt = "n")
      #   axis(1, at = c(1:M), labels = c(1:M), cex.axis = 1.2)
      #
      #   matlines(t(qn.dat[ip,]),type="b",col=c(4),ylim = c(low,up),xaxt = "n")
      #   matlines(t(mbqn.dat[ip,]),type="b",col=c(2),ylim = c(low,up),xaxt = "n")
      #   par(cex.lab=cex.lab) # is for y-axis
      #   par(cex.axis=cex.lab) # is for x-axis
      #
      # }else{
      #   plot(qn.dat[ip,],type="b",col=c(4),ylim = c(low,up),
      #        xlab = "sample",
      #        ylab = "normalized data", xaxt = "n",
      #        cex.lab = cex.lab,
      #        lwd  =1.5)
      #   axis(1, at = c(1:M), labels = c(1:M), cex.axis = 1)
      #   boxplot(mbqn.dat,col=(c("gold")),add = TRUE,notch=FALSE,
      #           main = "Median-balanced quantile normalized data\n & maximum RI/NRI features",
      #           cex.main = cex.main, outcex=0.3,xaxt = "n")
      #   lines(mbqn.dat[ip,],col=c(2),ylim = c(low,up), type="b", lwd=1.5)
      #
      # }
      # op <- par(cex = 0.8)
      # legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),
      #        col=c(4,2), lty=1, lwd = 1.5, cex=cex.legend, bty = "n") #fill = c(4,2),bty= "n",cex=1,ncol=1)

      if(save_fig){
        dev.copy2pdf(file=file.path(current.dir,fig4.name),width=8,height=4,out.type = "pdf")
        if(verbose) print(paste("Save figure to ",fig4.name))
      }
      }
      ###########################################################################################

      # # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn
      # # in one Figure
      # plot.new()
      # frame()
      #
      # if(is.null(filename)) {
      #   fig3.name <- "Figure_qn_nri_check.pdf"
      # }else{fig3.name <- paste0("Figure_qn_nri_check_", filename ,".pdf")}
      #
      # par(mfrow=c(2,1), mar = c(4,4,3,2))
      #
      # low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
      # up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))
      #
      # if(!is.null(feature_index)) ip <- unique(c(ip,feature_index))
      #
      # # if(length(ip)>1){
      #   boxplot(qn.dat,col=(c("gold")),notch=FALSE, ylim = c(low,up),
      #           xlab = "sample",
      #           ylab = "normalized data",
      #           main = "Quantile normalized data and maximum RI/NRI feature",
      #           cex.main = cex.main, outcex=0.4, xaxt = "n")
      #   axis(1, at = c(1:M), labels = c(1:M), cex.axis = .9)
      #   matlines(t(qn.dat[ip,]),type="b",pch=1, col=c(4),ylim = c(low,up), xaxt = "n")
      #   matlines(t(mbqn.dat[ip,]),type="b",pch=1,col=c(2),ylim = c(low,up), xaxt = "n")
      # }else{
      #   plot(qn.dat[ip,],col=c(4), type="b", ylim = c(low,up),xlab = "sample",
      #        ylab = "normalized data", xaxt = "n")
      #   axis(1, at = c(1:M), labels = c(1:M), cex.axis = .9)
      #   boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=FALSE,
      #           xlab = "sample", main = "Quantile normalized data and \n &maximum RI/NRI feature",
      #           cex.main = cex.main, outcex=0.4,xaxt = "n")
      #   lines(mbqn.dat[ip,],col=c(2),ylim = c(low,up), type="b", lwd=1.5, xaxt = "n")
      #
      # }
      # par(cex.lab=.8) # y-axis
      # par(cex.axis=.8) # x-axis
      # #op <- par(cex = 0.8)
      # legend(x = "bottomright",legend=(c("qn feature","mbqn feature")),col=c(4,2),
      #        lty=1, cex=0.8, bty = "n")
      # fig_label("A.", cex=1.7)

      # # ri/nri check and zero variance qn features
      # plot(as.table(res$p[1,])*100, xlim = xlimit, #xlim = c(0,N),
      #      xaxt = "n", col = c(1),
      #      ylim = c(0,100), xlab = "feature index",
      #      ylab = "frequency [%]",cex.main = .85,
      #      cex.lab = cex.lab, cex.axis = 1,cex.main =cex.main,
      #      main = "Occupation frequency of RI/NRI features\n & sample coverage of QN with zero variance")
      #
      # lines(not_nas*100,xlim = c(0,N),ylim = c(0,100), lty=1, col = c(2))
      #
      # xticklabels <- c(names(not_nas),colnames(p))
      # axis(1,at = as.integer(xticklabels), labels = FALSE, cex.axis = .8)
      # text(as.integer(xticklabels), par("usr")[3] - 5, labels = xticklabels,
      #      cex=0.75, srt = 90, pos = 1, xpd = TRUE)
      #
      # abline(h = 0, v = NA, col = "gray70")
      # abline(h = 100, v = NA, col = "gray70")
      # abline(h = 50, v = NA, col = "gray70")#,lty = "dashed")
      # abline(h = low_thr*100, v = NA, col = c(4),lty = "dashed")
      # op <- par(cex = .8)
      # legend(x = "topleft",inset = c(0,0.01),legend=c("RI/NRI features","coverage of RI features","threshold"),
      #        col=c(1,2,4), lty=c(1,1,2), cex=0.8, bty = "n")
      #
      # fig_label("B.", cex=2)

      # if(save_fig){
      #   dev.copy2pdf(file=file.path(current.dir,fig3.name),width=8,height=7,out.type = "pdf")
      #   if(verbose) print(paste("Save figure to ",fig3.name))
      # }
    } else { print("No NRI/RI present! You might want to adjust low_thr!")}

    } # fi show_fig

  return(invisible(res))

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

