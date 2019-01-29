#' Combined boxplot and line plot
#'
#' @param mtx a data matrix where rows represent features, e.g. protein abundance; columns represent groups, samples, e.g., replicates, experimental conditions.
#' @param irow an integer for row index or array of row indices for highlighting
#' @param vals numeric, array, matrix, or data frame. Array, matrix rows or data frame fields of length
#' \code{dim(mtx)[2]} are plot on top of the boxplot
# #' @param xlab xaxis-label
# #' @param ylab yaxis-label
# #' @param main figure title
#' @param filename save plot as pdf with filename in working directory
#' @param type line style for plot of row intensities
#' @param add.leg add legend to plot; default TRUE
#' @param ... additional arguments passed to plot, e.g. xlab, ylab, main, ylim, las
#' @details This function calls the \code{graphics::boxplot} function. Each box represents a group/column of the data matrix. Selected rows/features or user-defined arrays are plot as lines or points
#' on top of the boxes. Missing values are ignored.
#' @return Figure
#' @keywords quantile normalization, proteomics
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples ## Create boxplot of quantile normalized data matrix and plot
#' ## feature from median balanced quantile normalization on top of it.
#' \dontrun{
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3) # Create data matrix
#' qn.dat <- mbqn(x=X,FUN = NULL ,na.rm = TRUE) # Quantile normalization
#' mbqn.dat <- mbqn(x=X,FUN = median ,na.rm = TRUE) # Median balanced quantile normalization
#' ## Create boxplot and save output to file:
#' mbqnBoxplot(qn.dat,irow = 1, vals = mbqn.dat[1,], type = "b",filename = "fig_boxplot_qn.data.pdf")
#' }
#' @description Create a boxplot of a data matrix with groups. Plot selected features and/or additional user-defined data on top of it.
#' @importFrom grDevices dev.copy2pdf pdf
#' @importFrom graphics abline axis boxplot frame grconvertX grconvertY legend lines matlines par plot plot.new strheight strwidth text
#' @concept quantile, quantile normalization, rank invariance
#' @family data
#' @author Ariane Schad
#  August 2017
#' @export mbqnBoxplot
mbqnBoxplot <- function(mtx, irow = NULL, vals = NULL,filename = NULL, type = "l", add.leg = TRUE,...){

  if(!(is.matrix(mtx)|| is.data.frame(mtx))) {stop("Argument mtx must be a matrix or data.frame!")}

  opt.args <- list(...)
  xlab <- ifelse(is.null(opt.args$xlab), "sample", opt.args$xlab)
  ylab <- ifelse(is.null(opt.args$ylab), "intensity", opt.args$ylab)
  main <- ifelse(is.null(opt.args$main), "Boxplot",opt.args$main)
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

  #if(add.leg){
  # par(mar=c(par('mar')[1:3], 0)) # removes extraneous right inner margin space
  #}

  leg_text <- "data"
  lcol <- c("gold")
  lty <- 1
  las <- ifelse(length(opt.args$las)!=0, opt.args$las,0)


  if(length(irow)==1){
    leg_text <- c(leg_text,paste("protein",irow))
    lcol <- c("gold",2)
    lty <- c(lty,1)
  }else if(length(irow) >1){
    leg_text <- c(leg_text,paste("protein",irow))
    lcol <- c("gold",rep(2,length(irow)))
    lty <- c(lty,1:length(irow))
  }

  if(!is.null(vals)){

    if(class(vals)=="numeric" || class(vals) == "array"){
      lcol <- c(lcol,4)
      lty <- c(lty,1)
      leg.txt <- "feature"
    }

    if(class(vals) == "matrix" || class(vals) == "data.frame"){
      if(length(attributes(vals)$names)>=1){
        leg.txt <- as.array(names(vals))
      }else{
        leg.txt <- paste("feature",c(1:dim(vals)[1]))
      }
      lcol <- c(lcol,rep((1:dim(vals)[2])+3,each = 2)[1:min(dim(vals))])# seq(4,8)[1:min(dim(vals))])
      lty <- c(lty, rep(1:2, dim(vals)[2])[1:min(dim(vals))]) #c(lty,1:length(irow)) #rep(1,min(dim(vals)))) #seq(1:5)[1:min(dim(vals))])
    }
    leg_text <- c(leg_text,leg.txt)
  }

  if(add.leg){
    l <- legend(0, 0, bty='n', leg_text,
                plot=FALSE, pch=c(1, 2), lty=c(1, 2))
    # calculate right margin width in ndc
    w <- max(0.05,grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc'))
    par(omd=c(0, 1-w, 0, 1))
  }
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
  }else if(length(irow) >1){
    boxplot(mtx,col=(c("gold")),notch=F, plot = TRUE,
            xlab = xlab, ylab = ylab,
            main = main, ylim = ylim,
            las = las,...)
    matlines(t(mtx[irow,]),col=c(2),type = type, pch = 1,...)
  }

  if(!is.null(vals)){
    if(class(vals)=="array" || class(vals)=="numeric"){
      lines(vals,type=type, pch = 1,col=4,ylim = ylim,xlab = xlab, ylab = ylab,...)
    }else if(class(vals)=="matrix"){
      matlines(t(vals),type=type, pch = 1,lty = rep(1:2, dim(vals)[2]), col=rep((1:dim(vals)[2])+3,each = 2),...)
    }else{
      matlines(vals,type=type, pch = 1,lty = rep(1:2, dim(vals)[2]), col=rep((1:dim(vals)[2])+3,each = 2),...)
    }
  }

  if(add.leg){
    if(is.null(irow) & is.null(vals)){
      legend(par('usr')[2], par('usr')[4],
             bty='n', xpd=NA,
             legend=leg_text,
             fill = lcol[1],
             border = "black",
             col =  lcol,
             cex = cex.leg,
             box.lty=0,
             pt.cex = pt.cex)
    }else if(is.null(irow)){
      legend(par('usr')[2], par('usr')[4],
             bty='n', xpd=NA,
             legend=leg_text,
             fill = c(lcol[1],rep(NA,(length(lcol)-1))),
             border = c("black",rep(NA,(length(lcol)-1))),
             pch = rep(NA,(length(lcol))),
             lty = c(NA, lty[-1]),
             col =  lcol,
             cex = cex.leg,
             box.lty=0,
             pt.cex = pt.cex, y.intersp= y.intersp)
    }else{
      legend(par('usr')[2], par('usr')[4],
             bty='n', xpd=NA,
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
     print(paste("Save figure to",filename))
  }
}



