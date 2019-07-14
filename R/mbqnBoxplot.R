#' Combined box plot and line plot
#'
#' @description Create a box-and-whisker plot of a data matrix and
#' plot selected features and/or additional user-defined data on top of it.
#' @param mtx a matrix or data frame.
#' @param irow index or vector of row indices of matrix features to plot on top
#' of the boxplot.
#' @param vals numeric, array, matrix, or data frame of features with length
#' \code{ncol(mtx)} to plot on top of the boxplot.
#' @param add.leg add legend to plot.
#' @param ... additional arguments passed to the plot functions, e.g. xlab,
#' ylab, main, ylim, type, las.
#' @details This function calls \code{graphics::boxplot}.
#' Groups are represent by matrix columns. Selected rows/features or
#' user-defined arrays are plot on top of the box plot. Missing values are
#' ignored.
#' @return Figure.
#' @references Schad, A. and Kreutz, C., MBQN: R package for mean/median-
#' balanced quantile normalization. In prep. 2019
#' @examples ## Create boxplot of quantile normalized data matrix and plot
#' ## feature from median balanced quantile normalization on top of it.
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3) # Create data matrix
#' qn.dat <- mbqn(x=X,FUN = NULL ,na.rm = TRUE) # Quantile normalization
#' mbqn.dat <- mbqn(x=X,FUN = median ,na.rm = TRUE) # Median balanced quantile
#' normalization
#' ## Create boxplot:
#' plot.new()
#' mbqnBoxplot(qn.dat,irow = 1, vals = mbqn.dat[1,], type = "b")
#' @importFrom graphics axis boxplot grconvertX legend lines matlines par
#' strwidth
#' @family example
#' @author Ariane Schad
#  August 2017
#' @export mbqnBoxplot
mbqnBoxplot <- function(mtx,irow=NULL,vals=NULL,add.leg=TRUE, ...){
    filename = NULL
    if (!(is.matrix(mtx)|| is.data.frame(mtx))) {
        stop("Argument mtx must be a matrix or data.frame!")
    }
    opt.args <- list(...)
    type <- if (is.null(opt.args$type)) "l"
    cex.axis <- if (is.null(opt.args$cex.axis)) 0.8

    xlab <- ifelse(is.null(opt.args$xlab), "sample", opt.args$xlab)
    ylab <- ifelse(is.null(opt.args$ylab), "intensity", opt.args$ylab)
    main <- ifelse(is.null(opt.args$main), "Boxplot",opt.args$main)
    cex.leg <- ifelse(is.null(opt.args$cex), 0.8, opt.args$cex)
    cex <- ifelse(is.null(opt.args$cex), 0.8, opt.args$cex)
    pt.cex <- ifelse(is.null(opt.args$pt.cex), 0.8, opt.args$pt.cex)
    y.intersp <- ifelse(is.null(opt.args$y.intersp), 1, opt.args$y.intersp)
    fig.paper <- ifelse(is.null(opt.args$paper), "a4r", opt.args$paper)
    fig.width <- ifelse(is.null(opt.args$width), 10, opt.args$width)
    fig.height <- ifelse(is.null(opt.args$height), 5, opt.args$height)

    # y-axis range
    if (is.null(opt.args$ylim)){
        ylim <- range(mtx,na.rm = TRUE)
        if (!is.null(vals)){
            ymax <- max(ceiling(c(ylim,range(vals, na.rm = TRUE))))
            ymin <- min(floor(c(ylim,range(vals, na.rm = TRUE))))
        } else {
            ymax <- max(ceiling(ylim))
            ymin <- min(floor(ylim))
        }
        ymax <- ymax + 0.2*ymax
        ylim = c(ymin,ymax)
    } else {
        ylim <- opt.args$ylim
    }

    leg_text <- "data"
    lcol <- c("gold")
    lty <- 1
    las <- ifelse(length(opt.args$las)!=0, opt.args$las,0)

    if (length(irow)==1){
        leg_text <- c(leg_text,paste("id",irow))
        lcol <- c("gold",2)
        lty <- c(lty,1)
    } else if (length(irow) >1){
        leg_text <- c(leg_text,paste("id",irow))
        lcol <- c("gold",rep(2,length(irow)))
        lty <- c(lty,seq_len(length(irow)))
    }

    if (!is.null(vals)){
        if (is.numeric(vals) || is.array(vals)){
            lcol <- c(lcol,3)
            lty <- c(lty,1)
            leg.txt <- "feature"
        }
        if (is.matrix(vals) || is.data.frame(vals)){
            if (length(attributes(vals)$names)>=1){
                leg.txt <- as.array(names(vals))
            } else {
                leg.txt <- paste("feature",seq_len(nrow(vals)))
            }
            lcol <- c(
                lcol,rep(seq_len(ncol(vals))+2,each=6)[seq_len(min(dim(vals)))])
            lty <- c(lty, rep(seq_len(6), ncol(vals))[seq_len(min(dim(vals)))])
        }
        leg_text <- c(leg_text,leg.txt)
    }

    dy <- 0
    if (!is.null(colnames(mtx))) dy <- strwidth(colnames(mtx)[1],
                                                units = "figure",
                                                cex = cex)

    if(add.leg){
        #axis(side = 2, at = seq_len(18),labels = colnames(mtx), las =2)
        l <- legend(0, 0, bty='n', leg_text,
                    plot=FALSE, pch=c(1, 2), lty=c(1, 2), cex = cex.leg,
                    pt.cex = pt.cex,
                    y.intersp= y.intersp)
        # calculate right margin width in ndc
        w <- max(0.05,grconvertX(l$rect$w, to='ndc') - grconvertX(0, to='ndc'))
        par(omd=c(0, 1-w*.9, dy*3/4, 1))
    }

    if(!is.null(opt.args$ylim)) {
        opt.args <- .optargsReplace(..., replace = list(ylim = ylim))
    }

    if(!is.null(opt.args$width) ||
        !is.null(opt.args$height) ||
        !is.null(opt.args$y.intersp)) {
        opt.args <- .optargsRemove(..., remove=c("width","height","y.intersp"))
    }

    if (is.null(irow)){
        do.call(boxplot, c(list(x = mtx,use.cols = TRUE, col=c("gold"),
                                ylab = ylab,
                                xlab = xlab,
                                main = main,
                                cex = cex,
                                xlim = c(0,ncol(mtx)+.5),
                                las = las),opt.args))
    } else if (length(irow)==1){
        do.call(boxplot, c(list(x = mtx,use.cols = TRUE, col=c("gold"),
                                ylab = ylab,
                                xlab = xlab,
                                notch=FALSE,
                                main = main,
                                cex = cex,
                                xlim = c(0,ncol(mtx)+.5),
                                las = las),opt.args))
        do.call(lines, c(list(x=mtx[irow,], pch = 1,col=c(2)),opt.args))
    } else if (length(irow)>1){
        do.call(boxplot, c(list(x = mtx, use.cols = TRUE,
                                col=(c("gold")),
                                notch=FALSE, plot = TRUE,
                                xlab = xlab,
                                ylab = ylab,
                                main = main,
                                cex = cex,
                                xlim = c(0,ncol(mtx)+.5),
                                las = las),opt.args))
        do.call(matlines, c(list(y=t(mtx[irow,]),col=c(2), pch=1),opt.args))
    }

    if (!is.null(vals)){
        if (is.array(vals) || is.numeric(vals)){
            do.call(lines, c(list(x = vals, pch = 1,col=3,ylim=ylim,xlab=xlab,
                                ylab = ylab),opt.args))
        } else if (is.matrix(vals)){
            do.call(matlines, c(list(y = t(vals), pch = 1,
                lty = rep(seq_len(6), ncol(vals)),
                col=rep(seq_len(ncol(vals))+2,each = 6)),
                opt.args))
        } else { # data.frame
            do.call(matlines,
                c(list(y = vals,
                    pch = 1,
                    lty = rep(seq_len(6), ncol(vals)),
                    col=rep(seq_len(ncol(vals))+2,each = 6)),
                    opt.args))
        }
    }

    if (add.leg){
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
        } else if (is.null(irow)){
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
        } else {
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
}
