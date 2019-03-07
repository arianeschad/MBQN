#' Plot frequency of detected RI/NRI features
#'
#' @description Plot rank invariance frequency and feature coverage of detected RI and NRI features
#' @param obj list object of RI frequencies from \code{mbqnGetNRIfeatures()}.
#' @param save_fig logical indicating to save figure to pdf.
#' @param filename string to specify the filename of the .pdf.
#' @param verbose logical indicating to run function quietly.
#' @param ... additional arguments (cex, cex.lab, cex.axis, cex.main) passed to the plot function,
#' and (width, height, paper format(default "a4") of pdf) to \code{dev.copy2pdf()}.
#' @importFrom grDevices dev.copy2pdf
#' @details Graphical output of the NRI/RI identification results from \code{mbqnGetNRIfeatures()}.
#' For each detected NRI/RI feature, plot the feature index against the RI frequencies
#' together with the RI frequency detection threshold and print the sample coverage.
#' @return Figure, optionally saved in file Fig_nri_check_" \code{filename} ".pdf.
#' @seealso [mbqnGetNRIfeatures()] for detection of NRI/RI features.
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean/median-balanced quantile normalization. In prep. 2019
#' @examples ## Check data matrix for RI and NRI features
#' \dontrun{
#' x <- mbqnSimuData("omics.dep")
#' RI <- mbqnGetNRIfeatures(x, low_thr = 0.5, verbose = FALSE)
#' mbqnPlotRI(RI)
#' }
#' @author Ariane Schad
# 2017
#' @export mbqnPlotRI
mbqnPlotRI <- function(obj , save_fig = FALSE, filename = NULL,
                       verbose = FALSE, ...){

  if(!is.null(obj$nri)){

    ####### Graphical output #########

    # plot options
    opt.args <- list(...)
    cex.main  <- ifelse(is.null(opt.args$cex.main), 1.2, opt.args$cex.main)
    cex.legend <- ifelse(is.null(opt.args$cex), 0.8, opt.args$cex)
    cex.lab <- ifelse(is.null(opt.args$cex.lab), 1., opt.args$cex.lab)
    cex.axis <- ifelse(is.null(opt.args$cex.axis), 1., opt.args$cex.axis)
    y.intersp <- ifelse(is.null(opt.args$y.intersp), .8, opt.args$y.intersp)

    fig.paper <- ifelse(is.null(opt.args$paper), "a4", opt.args$paper)
    fig.width <- ifelse(is.null(opt.args$width), 10, opt.args$width)
    fig.height <- ifelse(is.null(opt.args$height), 5, opt.args$height)


    current.dir = getwd()

    par(mfrow=c(1,1), mar = c(4,4,3,2))

    if(is.null(filename)) {
      fig1.name <- "Figure_nri_check.pdf"
    }else{fig1.name <- paste0("fig_nri_check_", filename ,".pdf")}

    ylim <- c(min(as.numeric(names(obj$nri)))-1,max(as.numeric(names(obj$nri)))+2)

    dummy <- data.frame(frequency <- as.numeric(obj$p[1,]),
                        feature <- as.integer(colnames(obj$p)))
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
    abline(v = obj$low_thr*100, h = NA, col = c(4),lty = "dashed")
    axis(side = 2, at = dummy$`feature index`,cex.axis = 0.8)
    legend.txt <- cbind(legend.txt, "threshold")
    ind <- NULL


    coverage <- obj$p["sample.coverage",]

    par(mar = c(6,0.1,4,0))
    plot(NA, type="n", xlim=c(0,100),
         ylim=ylim,
         bty='n',
         xlab="",
         ylab="",
         yaxt="n",
         xaxt="n")
    if(length(coverage)>0){
      text(rep(5,length(coverage)), dummy$`feature index`,
           labels = paste0(as.character(coverage*100),"%"),
           col="red" ,cex = 0.8)
      legend.txt <- cbind(legend.txt, "feature coverage")
    }

    legend(x = 10, y = ylim[2],
           xpd = T,
           legend = legend.txt,
           col=c(1,4,2), lty=c(1,2,1), cex=cex.legend ,bty = "n", y.intersp = y.intersp)
    if(save_fig){
      dev.copy2pdf(file=file.path(current.dir,fig1.name),width=6,height=11,out.type = "pdf", paper="a4")
      if(verbose) print(paste("Save figure to ",fig1.name))
    }
  } else { print("No NRI/RI present! You might want to adjust low_thr!")}
}
