#' Plot RI/NRI feature frequencies and normalized/unnormalized features
#'
#' @description Check data matrix for rank invariant (RI) and
#' nearly rank invariant (NRI) features/rows across samples and visualize
#' result for different normalizations.
#' @param show_fig logical indicating whether results are displayed in figures.
#' @param save_fig logical indicating to save figures to pdf.
#' @param show_nri_only logical indicating to display and save only the RI/NRI detection graph to pdf.
#' @param filename a string for the .pdf-filenames.
#' @param ... additional plot arguments passed to \code{mbqnBoxplot}, \code{mbqnPlotRI} and \code{dev.copy2pdf()}.
#' @inheritParams mbqnGetNRIfeatures
#' @inheritParams mbqn
#' @inheritParams mbqnBoxplot
#' @inheritParams mbqnPlotRI
#' @importFrom grDevices dev.copy2pdf dev.off dev.size
#' @importFrom graphics abline layout plot text
#' @importFrom stats median
#' @details Rank data and check if lower and upper intensity tails are
#' dominated by few features. Apply quantile
#' normalization without and with mean-balancing and check the standard
#' deviation of normalized features located in the tails.
#' @return A set of figures that display the detected RI/NRI features and a list with elements:
#' \item{\code{p}}{a matrix with the rank invariance frequencies \code{ri.freq} and the sample coverage \code{sample.coverage} for all detected RI/NRI features}
#' \item{\code{max_p}}{maximum rank invariance frequency in percent}
#' \item{\code{ip}}{index of feature with maximum rank invariance frequency}
#' \item{\code{nri}}{table of the rank invariance frequencies in percent for each NRI/RI feature}
#' \item{\code{var0_feature}}{indices of features with zero sample variance after QN.}
#' @seealso [mbqnPlotRI()] and [mbqnBoxplot()] for the generation of figures, and [mbqn()] for normalization
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019
#' @examples ## Check data matrix for RI and NRI features
#' X <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA,1,3,1),ncol=3)
#' mbqnPlotAll(X, mean, low_thr = 0.5, save_fig = FALSE)
#' @author Ariane Schad
# 2017
#' @export mbqnPlotAll
mbqnPlotAll <- function(x, FUN = NULL,
                        low_thr = 0.5,
                        show_fig = TRUE,
                        save_fig = TRUE,
                        show_nri_only = FALSE,
                        filename = NULL,verbose = TRUE,...){

  opt.args <- list(...)

  res  <- mbqnGetNRIfeatures(x, FUN = FUN,
                             low_thr = low_thr,
                             verbose = verbose)

  # quantile normalisation and its standard deviation
  mbqn.dat <- mbqn(x = x, FUN = median, verbose = FALSE)
  qn.dat <- mbqn(x = x, FUN = NULL, verbose = FALSE)

  mbqn_ri.dat <- mbqnNRI(x = x, FUN = median, low_thr = low_thr, verbose = FALSE)

  ####### Graphical output #########

  if(show_fig){

    current.dir = getwd()

    # Occupation or rank invariance frequencies and sample coverage of RI/NRI features
    mbqnPlotRI(res, save_fig = save_fig, filename = filename,verbose = verbose, ...)

    # boxplot of quantile normalized data and maximum RI/NRI feature after qn and mbqn
    if(!show_nri_only){

      low <- floor(min(range(mbqn.dat,na.rm = TRUE)))
      up <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))

      df <- t(rbind(qn.dat[res$ip,],mbqn.dat[res$ip,]))
      colnames(df) <- c(paste0("QN",res$ip),paste0("MBQN",res$ip))
      df <- as.data.frame(df)

      fig2.name <- NULL
      if(save_fig){
        if(is.null(filename)) {
          fig2.name <- "Fig_example_qn_mbqn.pdf"
        }else{fig2.name <- paste0("Fig_example_qn_mbqn_", filename ,".pdf")}
        if(verbose) print(paste("Save figure to",fig2.name))
      }

      low2 <- floor(min(range(mbqn.dat,na.rm = TRUE)))
      up2 <- ceiling(max(range(mbqn.dat,na.rm = TRUE)))

      df2 <- t(rbind(qn.dat[res$ip,],mbqn.dat[res$ip,]))
      colnames(df2) <- c(paste0("QN",res$ip),paste0("MBQN",res$ip))
      df2 <- as.data.frame(df2)

      dev.off()
      plot.new()
      frame()
      mtx <- matrix(c(2, 1, 1, 3), byrow=TRUE, nrow=1)
      nf <- layout(mtx, heights=c(1), widths=c(6,5,1,0.5))
      # layout.show(nf)

      ylim2 <- c(low2,up2)
      ylim <- c(low,up)
      if(!is.null(opt.args$ylim)) {ylim <- ylim2 <- opt.args$ylim}


      opt.args.var <- opt.args

      if(!is.null(opt.args$ylim)) opt.args.var <- .optargsReplace(..., replace = list(ylim = ylim))
      opt.args.var <- .optargsRemove(opt.args.var, remove = c("main", "ylab"))


      # boxplot of mbqn data and with balanced qn RI/NRI features
      opt.args.var$ylab <- ""
      # remove empty list elements
      opt.args.var[which(lapply(opt.args.var, length)<1)] <- NULL
      do.call(mbqnBoxplot, c(list(mtx = mbqn.dat,
                                  vals = df2,
                                  add.leg = T,
                                  main = "MBQN data with RI/NRI feature"),
                             opt.args.var))

      # boxplot of qn data and with balanced qn RI/NRI features
      opt.args.var$ylab <- "normalized intensity"
      do.call(mbqnBoxplot, c(list(mtx = mbqn_ri.dat,
                                  vals = df,
                                  main = "QN data with RI/NRI feature",
                                  filename = fig2.name,
                                  add.leg = F),opt.args.var))


      ###########################################################################################

      # boxplot of qn data and with balanced qn RI/NRI features
      #fig4.name <- NULL
      #if(save_fig){
      #  if(is.null(filename)) {
      #    fig4.name <- "Figure_example_mbqn.pdf"
      #  }else{fig4.name <- paste0("Figure_example_mbqn_", filename ,".pdf")}
      #if(verbose) print(paste("Save figure to ",fig4.name))
      #}

    }
  }

  return(invisible(res))

}
