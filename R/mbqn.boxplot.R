#' Combined boxplot and line plot
#'
#' @param mtx data matrix where rows represent features, e.g. protein abundance; columns represent groups, samples, e.g., replicates, experimental conditions.
#' @param irow integer for row index or array of row indices for highlighting
#' @param vals numeric array/matrix of length \code{dim(mtx)[2]} that is plot on top of the boxplot
# #' @param xlab xaxis-label
# #' @param ylab yaxis-label
# #' @param main figure title
#' @param filename save plot to file with filename in working directory
#' @param type line style for plot of row intensities
# @inheritParams mbqn
#' @param add.leg add legend to plot; default TRUE
#' @param ... basic plot options like: xlab, ylab, main, ylim, las
#' @details This function calls the \code{graphics::boxplot} function. Each box represents a group/column of the data matrix. Selected rows/features or user-defined arrays are plot as lines or points
#' on top of the boxes. Missing values are ignored.
#' @return Figure
#' @keywords quantile normalization, proteomics
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
#' @examples ## Create boxplot of quantile normalized data matrix and plot
#' ## feature from median balanced quantile normalization on top of it.
# X <- matrix(rexp(20000, rate=.1), ncol=10)
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3) # Create data matrix
#' qn.dat <- mbqn(x=X,FUN = NULL ,na.rm = TRUE) # Quantile normalization
#' mbqn.dat <- mbqn(x=X,FUN = median ,na.rm = TRUE) # Median balanced quantile normalization
#' ## Create boxplot and save output to file:
#' mbqn.boxplot(qn.dat,irow = 1, vals = mbqn.dat[1,], type = "b",filename = "fig_boxplot_qn.data.pdf")
#' @description Create a boxplot of a data matrix with groups. Plot selected features and/or additional user-defined data on top of it.
#' @importFrom grDevices dev.copy2pdf dev.off dev.size pdf
#' @importFrom graphics abline axis boxplot frame grconvertX grconvertY legend lines matlines par plot plot.new strheight strwidth text
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @export mbqn.boxplot
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#  August 2017
#' @export mbqn.boxplot
mbqn.boxplot <- function(mtx, irow = NULL, vals = NULL,filename = NULL, type = "l", add.leg = TRUE,...){
  #xlab = NULL, ylab = NULL, main = NULL

  if(!is.matrix(mtx)) {warning("Data must be a matrix!"); break}

  opt.args <- list(...)
  xlab <- ifelse(is.null(opt.args$xlab), "sample", opt.args$xlab)
  ylab <- ifelse(is.null(opt.args$ylab), "intensity", opt.args$ylab)
  main <- ifelse(is.null(opt.args$main), "Boxplot data matrix",opt.args$main)
  cex.leg <- ifelse(is.null(opt.args$cex), 0.8, opt.args$cex)
  pt.cex <- ifelse(is.null(opt.args$pt.cex), 0.8, opt.args$pt.cex)
  cex.axis <- ifelse(is.null(opt.args$cex.axis), 0.8, opt.args$cex.axis)
  y.intersp <- ifelse(is.null(opt.args$y.intersp), 1, opt.args$y.intersp)

  if(!is.null(filename)){
    pdf(paste0(filename,".pdf"), width=10,height=6,paper="a4r")
  }

  # y-axis range
  if(is.null(opt.args$ylim)){
    ylim <- range(mtx,na.rm = TRUE)
    if(!is.null(vals)){
      ymax <- max(ceiling(c(ylim,range(vals, na.rm = T))))
      ymin <- min(floor(c(ylim,range(vals, na.rm = T))))
    }else{
      ymax <- max(ceiling(ylim))
      ymin <- min(floor(ylim))
    }
    ymax <- ymax + 0.2*ymax
    ylim = c(ymin,ymax)
  }else{
    ylim <- opt.args$ylim}


  # Set plot layout
 # if(max(par()$mfrow)==1) {
  layout(mat = matrix(c(1,2), nrow = 1, ncol = 2), widths = c(5, 2))# Widths of the two columns

  # set graphic window
  par(mar=c(6.9, 4.1, 4.1, .2), cex.axis = cex.axis, no.readonly = TRUE)

 # }
  #par(new = TRUE, mar=c(6.9, 4.1, 4.1, 9.1), xpd=TRUE, cex.axis = 0.8)
  # grid(nx=NA, ny=NULL) #grid over boxplot

  leg_text <- "data"
  lcol <- c("gold")
  lty <- 1
  las <- ifelse(length(opt.args$las)!=0, opt.args$las,0)

  if(is.null(irow)){
    boxplot(mtx,col=c("gold"),
            ylab = ylab, xlab = xlab,
            main = main, xlim = c(0,dim(mtx)[2]+.5),
            cex = 0.8,
            ylim = ylim,
            las = las,...)
  }else if(length(irow)==1){
    boxplot(mtx,col=(c("gold")),notch=F,
            xlab = xlab, ylab = ylab,
            main = main,
            ylim = ylim,
            las = las,...)
    lines(mtx[irow,],type = type, pch = 1,col=c(2),...)
    leg_text <- c(leg_text,paste("protein",irow))
    lcol <- c("gold",2)
    lty <- c(lty,1)
  }else if(length(irow) >1){
    boxplot(mtx,col=(c("gold")),notch=F, plot = TRUE,
            xlab = xlab, ylab = ylab,
            main = main, ylim = ylim,
            las = las,...)
    matlines(t(mtx[irow,]),col=c(2),type = type, pch = 1,...)
    leg_text <- c(leg_text,paste("protein",irow))
    lcol <- c("gold",rep(2,length(irow)))
    lty <- c(lty,1:length(irow))
  }

   if(!is.null(vals)){
     if(length(attributes(vals)$names)>=1){
       #leg.txt <- "feature"}else{
       leg.txt <- as.array(names(vals))
     }else{leg.txt <- "feature"}

  #  #leg.bxplt <- ifelse(length(names(vals))>0, names(vals), "feature")
  #  if(is.data.frame(vals)) vals <- as.matrix(vals)

  #  # leg.txt <- ifelse(is.null(colnames(vals)), ifelse(is.null(names(vals)), "feature", names(vals)), colnames(vals))
  #  if(is.null(colnames(vals))){
  #    if(is.null(colnames(vals))){
  #      leg.txt <- "feature"
  #    }else{
  #     leg.txt <- names(vals)
  #     }

if(is.null(dim(vals))){
  lines(vals,type=type, pch = 1,col=4,ylim = ylim,xlab = xlab, ylab = ylab,...)
  lcol <- c(lcol,4)
  lty <- c(lty,1)
}else{
  matlines(vals,type=type, pch = 1,lty = rep(1:2, dim(vals)[2]), col=rep((1:dim(vals)[2])+3,each = 2),...)
  lcol <- c(lcol,rep((1:dim(vals)[2])+3,each = 2)[1:min(dim(vals))])# seq(4,8)[1:min(dim(vals))])
  lty <- c(lty, rep(1:2, dim(vals)[2])[1:min(dim(vals))]) #c(lty,1:length(irow)) #rep(1,min(dim(vals)))) #seq(1:5)[1:min(dim(vals))])
  }
  leg_text <- c(leg_text,leg.txt)

   }

#par(mar =c(6.9, 0.1, 4.1, .2), cex.axis =0.8)#, mfrow = par()$mfrow)
#if(max(par()$mfrow)==1)
  par(mar=c(6.9, 0.1, 4.1, .2), cex.axis = cex.axis, no.readonly = TRUE)#, mfrow = par()$mfrow)

if(add.leg){
plot(NA, type="n",
     bty="n", yaxt="n",cex = 0.8,
     ylim = ylim,
     xaxt="n", xlab = '', ylab = '')
if(is.null(irow)){
  legend(x = "topleft",
         legend=leg_text,
         #fill = c(lcol[1],NA),
         #lty = lty[-1]),
         #lty = c(NA, lty[2]),
         #pch = c(NA,NA),
         #border = c("black",NA),
         fill = c(lcol[1],rep(NA,(length(lcol)-1))),
         border = c("black",rep(NA,(length(lcol)-1))),
         pch = rep(NA,(length(lcol))),
         lty = c(NA, lty[-1]),
         col =  lcol,
         cex = cex.leg,
         box.lty=0,
         pt.cex = pt.cex, y.intersp= y.intersp)
}else{
  legend(x = "topleft",
         # y.intersp = 0.5,  # vertical spacing between legend entries
         #xpd = TRUE,
         # legend=leg_text[-1],
         # lty = lty[-1],
         # col=lcol[-1],
         # cex = cex.leg,
         # box.lty=0,
         # pt.cex = pt.cex,
         # y.intersp= y.intersp
         legend=leg_text,
         fill = c(lcol[1],rep(NA,(length(lcol)-1))),
         border = c("black",rep(NA,(length(lcol)-1))),
         pch = rep(NA,(length(lcol))),
         lty = c(NA, lty[-1]),
         col=lcol,
         cex = cex.leg,
         box.lty=0,
         pt.cex = pt.cex,
         y.intersp= y.intersp
  )
}
}

if(!is.null(filename)){
  dev.off()
  print(paste("Save figure to",filename))
}
}



