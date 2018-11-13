#' Boxplot of data matrix with selected rows highlighted
#'
#' @param mtx A data matrix. Rows represent features, e.g. protein abundances; columns represent groups, samples, e.g., experimental conditions, replicates.
#' @param irow row index or array of row indices for highlighting
#' @param vals numeric array of length \code{dim(x)[2]} that is plot on top of boxplot
#' @param xlab xaxis-label
#' @param ylab yaxis-label
#' @param main figure title
#' @param filename save plot to file with filename in working directory
#' @param type line style for plot of row intensities
# @inheritParams mbqn
#' @details Create a boxplot of a data matrix and highlight values of selected rows with lines or points
#' across samples. Plot additional user-defined lines on top of the boxplot.
#' @return Figure
#' @keywords quantile normalization, proteomics
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
#' @examples ## Creat boxplot of quantile normalized data matrix, plot feature from median balanced quantile normalization on top of it.
#' X <- matrix(rexp(20000, rate=.1), ncol=10)
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3) # Create data matrix
#' qn.dat <- mbqn(x=X,FUN = NULL ,na.rm = TRUE) # Quantile normalization
#' mbqn.dat <- mbqn(x=X,FUN = median ,na.rm = TRUE) # Median balanced quantile normalization
#' mbqn.boxplot(qn.dat,irow = 1, vals = mbqn.dat[1,], type = "b",filename = "fig_boxplot_qn.data.pdf") # Create boxplot and save output to file
#' @description Create a boxplot of data matrix and highlight selected features. Plot additional user-defined data on top of it.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#  August 2017
#' @export
mbqn.boxplot <- function(mtx, irow = NULL, vals = NULL, xlab = NULL, ylab = NULL, main = NULL,filename = NULL, type = "l"){

  if(!is.matrix(mtx)) {warning("Data must be a matrix!"); break}


  if(is.null(xlab)) xlab <- "sample"
  if(is.null(ylab)) ylab <- "intensity"
  if(is.null(main)) main <- "Boxplot data matrix"

  if(!is.null(filename)){
    pdf(paste0(filename,".pdf"), width=6,height=5,paper='special')
  }

  # y-axis range
  #if(is.null(ylim)){
  ylim <- range(mtx,na.rm = T)
  if(!is.null(vals)){
    ymax <- max(ceiling(c(ylim,range(vals, na.rm = T))))
    ymin <- min(floor(c(ylim,range(vals, na.rm = T))))
  }else{
    ymax <- max(ceiling(ylim))
    ymin <- min(floor(ylim))
  }
  ymax <- ymax + 0.2*ymax
  ylim = c(ymin,ymax)
  #}

  # set a new graphic window
  frame()
  # Add extra space to right of plot area; change clipping to figure
  par(new = TRUE, mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE, cex.axis = 0.8)

  #grid(nx=NA, ny=NULL) #grid over boxplot

  leg_text <- "data"
  lcol <- c("gold")
  lty <- 1
  if(is.null(irow)){
    boxplot(mtx,col=c("gold"),ylab = ylab, xlab = xlab, main = main, xlim = c(0,dim(mtx)[2]+.5),cex = 0.8)
  } else if(length(irow)==1){
    boxplot(mtx,col=(c("gold")),notch=F, xlab = xlab, ylab = ylab, main = main)
    lines(mtx[irow,],type = type, pch = 1,col=c(2),ylim = ylim)
    leg_text <- c(leg_text,paste("protein",irow))
    lcol <- c("gold",2)
    lty <- c(lty,1)
  } else if(length(irow) >1){
    boxplot(mtx,col=(c("gold")),notch=F, plot = TRUE, xlab = xlab, ylab = ylab, main = main)
    matlines(t(mtx[irow,]),col=c(2),type = type, pch = 1, ylim = ylim)
    leg_text <- c(leg_text,paste("protein",irow))
    lcol <- c("gold",rep(2,length(irow)))
    lty <- c(lty,1:length(irow))
  }
  if(!is.null(vals)){
    n <- if(is.null(names(vals))) "array" else names(vals)
    if(is.data.frame(vals))  vals <- data.matrix(vals)
    lines(vals,type=type, pch = 1,col=c(4),ylim = ylim,xlab = xlab, ylab = ylab)
    leg_text <- c(leg_text,n)
    lcol <- c(lcol,4)
    lty <- c(lty,1)
  }

  if(is.null(irow)){
    legend(x = "topright", inset=c(-0.3,0),legend=leg_text, fill = lcol, col =  lcol, cex =.8, box.lty=0, pt.cex = 0.8)
  }else{
    legend(x = "topright", inset=c(-0.4,0),
           legend=leg_text[-1], lty = lty[-1],
           col=lcol[-1], cex =.8, box.lty=0, pt.cex = 0.8)
  }

  if(!is.null(filename)){
    dev.off()
    print(paste("Save figure to",filename))
  }

}
