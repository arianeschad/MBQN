#' Example illustrating normalization of LFQ protein itensities with MBQN
#'
#' @description This function illustrates the check for rank invariant (RI)
#' and nearly rank invariant (NRI) features and the normalization of LFQ intensities
#' of proteins for selected datasets deposited to PRIDE \[1\].
#' @param which.example numerical value between 1-4 to select one of four datasets.
#' @param source.path character specifying the path where to store and search for the
#' "proteinGroups.txt"-file; default = NULL uses current working directory.
#' @inheritParams mbqn
#' @importFrom grDevices dev.copy2pdf dev.off
#' @importFrom stats median
#' @return Figures saved as pdf-files.
#' @details To run this function, download the proteinGroups.txt file
#' of the selected example manually in advance or use the automatic dowload
#' provided by the function. \cr\cr
#' Datasets available from PRIDE supported by this function:\cr
#' 1 PXD001584 - contains one RI feature (default) (file size 10.8 MB) c\cr
#' 2 PXD005138 - contains one RI feature (file size 7440.6 MB) \[3\]\cr
#' 3 PXD005861 - contains one RI feature (file size 334.7 MB) \[4\]\cr
#' 4 PXD006617 - contains one RI feature (file size 290.8 MB) \[5\]\cr\cr
#' Check the reference for further details on the examples.
#' These examples are also listed in \[6\].
# @return A matrix of median- or mean-balanced quantile normalized data
#' @family example
#' @seealso [get_pxdfile()] for download dataset from PRIDE and [loadFile()] for reading LFQ intensities from file.
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I, Mayer G,
#' Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.
#' \[2\] Ramond E, et al. Importance of host cell arginine uptake in Francisella
#' phagosomal escape and ribosomal protein amounts. Mol Cell Proteomics. 2015 14(4):870-881\cr
#' \[3\] Beaumont V, et al. Phosphodiesterase 10A Inhibition Improves Cortico-Basal
#' Ganglia Function in Huntington's Disease Models. Neuron. 2016 92(6):1220-1237\cr
#' \[4\] Turetschek R, Desalegn G, Epple T, Kaul HP, Wienkoop S. Key metabolic
#' traits of Pisum sativum maintain cell vitality during Didymella pinodes infection:
#' cultivar resistance and the microsymbionts' influence. J Proteomics.
#' 2017 Mar 4. pii: S1874-3919(17)30075-1\cr
#' \[5\] Ranjbar Sistani N, Kaul HP, Desalegn G, Wienkoop S. Rhizobium Impacts
#' on Seed Productivity, Quality, and Protection of Pisum sativum upon Disease
#' Stress Caused by Didymella pinodes: Phenotypic, Proteomic, and Metabolomic
#' Traits. Front Plant Sci. 2017 8:1961\cr
#' \[6\] Schad, A. and Kreutz, C., MBQN: R package for mean/median-balanced quantile
#' normalization. In prep. 2019\cr
#' @examples ## Check LFQ intensities of the protein dataset in
#' ## PXD001584 for RI and NRI features
#'\dontrun{
#' mbqnExample(which.example = 3)
#'}
#' @author Ariane Schad
# 2017
#' @export mbqnExample
mbqnExample <- function(which.example = NULL, source.path = NULL){


  # PXD repositories
  ids <- c("PXD001584","PXD005138","PXD005861","PXD006617")

  if(is.null(which.example) | !(which.example %in% c(1,2,3,4))) stop("Error: Select an example between 1-4")

  pxd_id <- ids[which.example]

  # Load file
  out <- loadFile(pxd_id, source.path = source.path, file.pattern = "proteinGroups.txt")
  featureAnnotations <- out$featureAnnotations
  mtx <- out$mtx

  # filter for potential contaminants and identified only by site features
  mtx <- mtx[-out$ixs,]
  featureAnnotations <- featureAnnotations[-out$ixs,]

  low_thr <- 0.5
  ylim <- NULL
  ix <- seq(1,ncol(mtx))

  if(pxd_id == "PXD001584") {
    ix <- c(seq(1,9), seq(19,27))
    mtx <- mtx[,ix]
    ylim.qn <- ylim <- c(22.5,36)
  }

  if(pxd_id == "PXD006617") {
    ylim.qn <- ylim <- c(16.5,36)
    low_thr <- 0.5
  }

  ####################################################################################

  res <- mbqnPlotAll(mtx,
                     FUN = median,
                     low_thr = low_thr,
                     las = 2,
                     type = "l",
                     feature_index = NULL,
                     show_fig = TRUE,
                     filename = pxd_id,
                     show_nri_only = TRUE,
                     save_fig = TRUE,
                     axis.cex = 0.5,
                     y.intersp= 0.5)

  dev.off()

  # get protein name of strongest nri/ri feature
  nri_max <- as.numeric(names(which.max(res$nri)))
  #featureAnnotations$proteinDescription[nri_max]
  #featureAnnotations$proteinName[nri_max]
  #featureAnnotations$nbPeptides[nri_max]

  df <- lapply(seq(1,dim(featureAnnotations)[2]),function(j) featureAnnotations[[j]][[nri_max]])
  names(df)<- names(featureAnnotations)

  # truncate long name strings
  df$proteinName <-  paste(strtrim(df$proteinName,80),"...")
  df$proteinDescription <- paste(strtrim(df$proteinDescription,80),"...")

  colnames(mtx) <- gsub("LFQ intensity","",colnames(mtx))
  mbqn.mtx <- mbqn(mtx,FUN = median)
  qn.mtx <- mbqn(mtx,FUN = NULL)

  # Boxplot of QN intensity features, highlight RI/NRI Features
  save.fig <- TRUE
  fig1.name <- paste0("fig_qnLFQ_", pxd_id,".pdf")


  if(length(ylim)==0) {
    ylim <- c(floor(min(range(mbqn.mtx, na.rm = TRUE))),ceiling(max(range(mbqn.mtx, na.rm = TRUE))))
    ylim.qn <- c(floor(min(range(qn.mtx, na.rm = TRUE))),ceiling(max(range(mbqn.mtx, na.rm = TRUE))))
  }
  #colnames(qn.mtx) <- colnames(mbqn.mtx) <- ix

  plot.new()
  frame()
  mbqnBoxplot(mtx = qn.mtx,main = paste("QN"),
              ylab = "LFQ intensity",
              ylim = ylim, xlab = "",las =2,
              irow = c(as.numeric(names(res$nri))), y.intersp = 0.5)
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig1.name),width=10,height=6,paper="a4r",out.type = "pdf")
    print(paste("Save figure to ",fig1.name))
  }

  #if(which.example==1){
    fig3.name <- paste0("fig_mbqn_vs_qn_", pxd_id,".pdf")
    plot.new()
    frame()

    # select a qn feature:
    is.full.feature <- which((apply(is.na(mtx),1,sum))==0)[1]
    df <- data.frame(QN.feature = qn.mtx[is.full.feature,])
    names(df) <- paste("QN feature", is.full.feature)
   # colnames(mbqn.mtx ) <- colnames(mtx)
    mbqnBoxplot(mtx = mbqn.mtx,main = "MBQN with QN-feature",
                irow = c(as.numeric(names(res$nri)),is.full.feature),
                xlab= "", las =2,
                ylab = "LFQ intensity",
                vals = df,
                ylim = ylim, y.intersp = 0.5)
    if(save.fig){
      dev.copy2pdf(file=file.path(getwd(),fig3.name),width=10,height=6,paper="a4r", out.type = "pdf")
      print(paste("Save figure to ",fig3.name))
    }
  #}

  # QN of data and balance only NRI/RI features
  fig4.name <- paste0("fig_qn_balance_nri_only_", pxd_id,".pdf")
  dev.off()
  plot.new()
  frame()
  par(mfrow = c(1,1))


  df <- data.frame(qn.mtx[as.numeric(names(res$nri)),])
  if(ncol(df)==1) df <- t(df)
  rownames(df) <- paste("QN feature",names(res$nri))
  df2 <- data.frame(mtx[as.numeric(names(res$nri)),])
  if(ncol(df2)==1) df2 <- t(df2)
  rownames(df2) <- paste("unnormal. feature",names(res$nri))
  colnames(df) <- colnames(df2)
  df <- rbind(df,df2)
  df <- as.data.frame(t(df))

  mtx.nri <- mbqnNRI(mtx,FUN = median,low_thr = 0.5, verbose = FALSE)
  mbqnBoxplot(mtx = mtx.nri,
              irow = as.numeric(names(res$nri)),
              ylim = ylim.qn, xlab= "", las=2,
              ylab = "LFQ intensity",
              vals = df,lwd = 1.,
              main = "QN with RI/NRI balanced",
              cex.axis = 1, cex.lab = .9, cex = .9, y.intersp = 0.5)

  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig4.name),width=9,height=6,paper="a4r", out.type = "pdf")
    print(paste("Save figure to ",fig4.name))
  }
}
