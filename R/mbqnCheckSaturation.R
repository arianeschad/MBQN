#' Check data matrix for rank invariant (RI) and nearly rank invariant (NRI) features
#'
#' @param dat A data matrix. Rows - features, e.g. protein abundances; columns - samples
#' @param FUN median, mean, or another function used to balance features across colums. If left empty, quantile normalization
#' is applied without balancing the data
#' @param low_thr Numerical value for the lower threshold for NRI frequency, default = 0.5
#' @param feature_index Integer that indicates the index of a feature of interest that is plotted in the boxplot; default NULL
#' @param show_fig Logical flag indicating whether results should be displayed; default = TRUE.
#' @param save_fig Logical to save figures to file
#' @param show_nri_only Logical to print and save only the RI/NRI detection graph; default FALSE
#' @param filename String for naming figures, default = NULL
#' @param verbose Logical for running function quiet
#' @param ... Optional plot arguments passed to \code{mbqn.boxplot}
#' @inheritParams mbqn
#' @inheritParams mbqnBoxplot
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
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note., 2018
#' @examples ## Check data matrix for RI and NRI features
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqnCheckSaturation(X, mean, low_thr = 0.5, save_fig = FALSE)
#' @description Check data matrix for intensity features in rows which dominate the upper tail, i.e. for
#' features that have constant or nearly constant rank across samples.
#' @author A. Schad, \email{ariane.schad@zbsa.de}
# 2017
#' @export mbqnCheckSaturation
mbqnCheckSaturation <- function(dat, FUN = NULL,
                                low_thr = 0.5,
                                feature_index = NULL,
                                show_fig = TRUE,
                                save_fig = TRUE,
                                show_nri_only = FALSE,
                                filename = NULL,verbose = TRUE,...){


  res  <- mbqnGetNRIfeatures(dat, FUN = FUN,
                             low_thr = low_thr,
                             verbose = verbose)

  # quantile normalisation and its standard deviation
  qn.dat <- mbqn(x = dat,FUN = NULL, verbose = FALSE)
  mbqn.dat <- mbqnNRI(x = dat, FUN = median, low_thr = low_thr, verbose = FALSE)

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

        ind <- as.numeric(names(res$nri))
      if(length(ind)>0) {
        if(length(ind)==1) {
          coverage <- sum(!is.na(qn.dat[ind,]))/dim(qn.dat)[2]
        } else { coverage <- apply(!is.na(qn.dat[ind,]), 1,sum)/dim(qn.dat)[2]
        }
      }


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
             legend = legend.txt,
             col=c(1,4,2), lty=c(1,2,1), cex=cex.legend ,bty = "n", y.intersp = 0.3)
      if(save_fig){
        dev.copy2pdf(file=file.path(current.dir,fig1.name),width=6,height=11,out.type = "pdf", paper="a4")
        if(verbose) print(paste("Save figure to ",fig1.name))
      }

      ###########################################################################################
      # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn
      if(!show_nri_only){
        dev.off()
        plot.new()
        frame()
        #par(mfrow=c(1,1), mar = c(4,4,3,2))


        low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
        up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))

        mbqnBoxplot(qn.dat,
                    vals = data.frame(QN.feature = qn.dat[res$ip,], MBQN.feature = mbqn.dat[res$ip,]),
                    ylim = c(low,up),
                    ylab = "normalized intensity",
                    main = "QN data with unbalanced \n and balanced maximum RI/NRI feature",
                    cex.main = 0.8,...)

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

        mbqnBoxplot(mbqn.dat,
                    vals = data.frame(QN.feature = qn.dat[res$ip,], MBQN.feature = mbqn.dat[res$ip,]),
                    ylim = c(low,up),
                    ylab = "normalized intensity",
                    main = "MBQN data with unbalanced \n and balanced maximum RI/NRI feature",
                    cex.main = 0.8,...)

        if(is.null(filename)) {
          fig4.name <- "Figure_example_mbqn.pdf"
        }else{fig4.name <- paste0("Figure_example_mbqn_", filename ,".pdf")}

        if(save_fig){
          dev.copy2pdf(file=file.path(current.dir,fig4.name),width=8,height=4,out.type = "pdf")
          if(verbose) print(paste("Save figure to ",fig4.name))
        }
      }
    } else { print("No NRI/RI present! You might want to adjust low_thr!")}

  } # fi show_fig

  return(invisible(res))

}
