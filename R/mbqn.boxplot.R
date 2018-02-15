#' Boxplot of data matrix with highlighted selected row intensities
#'
#' @param x A data matrix. Rows represent features, e.g. protein abundances; columns represent groups, samples, e.g., experimental conditions, replicates.
# @param FUM mean or median, if left empty, quantile normalization
# is applied without balancing the data
#' @param irow row index or array of row indices for highlighting
#' @param vals numeric array of length \code{dim(x)[2]} to plot on top of boxplot
#' @param xlab xaxis-label
#' @param ylab yaxis-label
#' @param main figure title
#' @param filename save plot to file with filename
#' @param type line style for plot of row intensities
# @inheritParams mbqn
#' @details Create boxplot of data and highlight values of selected rows with lines or points
#' across samples
#' @return Figure
#' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C. (2017) Mean-balanced quantile
#' normalization for processing label-free quantitative proteomics
#' data with abundance-isolated proteins. Biostatistics xxx in prep.
#' @examples qn.dat <- mbqn(x=dat,FUN = NULL ,na.rm = TRUE)
#' mbqn.dat <- mbqn(x=dat,FUN = median ,na.rm = TRUE)
#' mbqn.boxplot(qn.dat,irow = 1, vals = mbqn.dat[1,], filename = "fig_boxplot_qn.data.pdf")
#' @description Boxplot of a data matrix dominate and highlight selected row indices.
#' @author A. Schad, \email{ariane.schad@zbsa.de}


mbqn.boxplot <- function(x, irow = NULL, vals = NULL, xlab = NULL, ylab = NULL, main = NULL,filename = NULL, type = "l"){

  if(!is.matrix(x)) {warning("Data must be a matrix!"); break}


  if(is.null(xlab)) xlab <- "sample"
  if(is.null(ylab)) ylab <- "data"
  #if(is.null(main)) main <- "Normalized data and maximum value"

  if(!is.null(filename)){
    pdf(paste0(filename,".pdf"), width=6,height=5,paper='special')
  }


  # y-axis range
  ylim <- range(x)
  ymax <- max(ceiling(c(ylim,range(vals))))
  ymax <- ymax + 0.2*ymax
  ymin <- min(floor(c(ylim,range(vals))))
  ylim = c(ymin,ymax)

  #ymax <- ceiling(max(max(qn.dat),max(mbqn.dat)))
  #ymax <- ymax + 0.2*ymax
  #ymin <- floor(min(min(qn.dat),min(mbqn.dat)))

  # if(length(irow)==1){
  #   plot(qn.dat[irow,],type=type,col=c(4),ylim = ylim,xlab = xlab, ylab = ylab)
  #   boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = xlab, main = main)
  #   lines(mbqn.dat[irow,],type=type,col=c(3),ylim = ylim,xlab = xlab, ylab = ylab)
  # } else if(length(irow) >1){
  #   matplot(t(qn.dat[irow,]),type = type,col=c(4),ylim = ylim, xlab = xlab,ylab = ylab)
  #   boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = xlab, main = main)
  #   matlines(t(mbqn.dat[irow,]),type=type,col=c(3),ylim = ylim, xlab = xlab,ylab = ylab)
  # } else if(is.null(irow)){
  #   boxplot(qn.dat,col=(c("gold")),add = TRUE,notch=TRUE, xlab = xlab, main = main)
  #   matlines(t(mbqn.dat[irow,]),type=type,col=c(3),ylim = ylim, xlab = xlab,ylab = ylab)
  # }

  # set a new graphic window
  frame()
  grid(nx=NA, ny=NULL) #grid over boxplot
  par(new=TRUE)

  leg_text <- "data"
  lcol <- c("gold")
  lty <- 1
  if(is.null(irow)){
    boxplot(x,col=c("gold"),ylab = ylab, xlab = xlab, main = main, xlim = c(0,dim(x)[2]+.5))
  } else if(length(irow)==1){
    boxplot(x,col=(c("gold")),notch=TRUE, xlab = xlab, main = main)
    lines(x[irow,],type= type ,col=c(3),ylim = ylim,xlab = xlab, ylab = ylab)
    leg_text <- c(leg_text,paste("row",irow))
    lcol <- c("gold",3)
    lty <- c(lty,1)
  } else if(length(irow) >1){
    boxplot(x,col=(c("gold")),notch=TRUE, plot = TRUE, xlab = xlab, main = main)
    matlines(t(x[irow,]),type=type,col=c(3),ylim = ylim, xlab = xlab,ylab = ylab)
    leg_text <- c(leg_text,paste("row",irow))
    lcol <- c("gold",rep(3,length(irow)))
    lty <- c(lty,1:length(irow))
  }
  if(!is.null(vals)){
    n <- if(is.null(names(vals))) "array" else names(vals)
    if(is.data.frame(vals))  vals <- data.matrix(vals)
    lines(vals,type=type,col=c(2),ylim = ylim,xlab = xlab, ylab = ylab)
    leg_text <- c(leg_text,n)
    lcol <- c(lcol,2)
    lty <- c(lty,1)
  }

  legend(x = "topright", inset=.02,legend=leg_text, lty = lty, col =  lcol, cex =1, box.lty=0)
  #legend(x = "topright",legend=(c("qn","mbqn")),fill = c(4,3),bty= "n",cex=1,ncol=1)

  if(!is.null(filename)){
    dev.off()
    print(paste("Save figure to",filename))
  }

}
