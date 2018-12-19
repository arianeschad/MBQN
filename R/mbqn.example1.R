#' Example illustrating the analysis of PRIDE data with the MBQN package
#'
#' @description This function runs an example to illustrate the check of LFQ intensities
#' of proteins for rank invariant (RI) and nearly rank invariant (NRI) features
#' and the normalization of the data. It downloads a selected dataset from PRIDE
#' or it uses the dataset in filepath. The example needs the R package "rpx" to
#' download data from PRIDE.
#' @param which.example Numerical values between 1-4, default = 1,
#' to select between 4 dataset.
#' @param source.path Character with path where to store and search for proteinGroups.txt
#' file; default = NULL uses "installationpath/MBQN/examples".
#' @inheritParams mbqn
#' @importFrom grDevices dev.copy2pdf dev.off dev.size pdf
#' @importFrom filesstrings file.move
#' @details Collecting information on the experiments requires the package rpx by Laurent Gatto (2017),
#' version 1.10.2 from https://github.com/lgatto/rpx.\cr
# "rpx: R Interface to the ProteomeXchange Repository".
# R package version
#' The function used to read the proteinGroups.txt files uses source code
#' of SafeQuant::parseMaxQuantProteinGroupTxt, \cr
# #' See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.1.
#' version 2.3.1, https://CRAN.R-project.org/package=SafeQuant.\cr\cr
#' Experiments on PRIDE that are supported by this Example:\cr
#' 1 PXD001584 - contains a RI (default) (file size 10.8 MB) \[1\]\cr
#' 2 PXD005138 - contains a RI feature (file size 7440.6 MB) \[2\]\cr
#' 3 PXD005861 - contains RI feature (file size 334.7 MB) \[3\]\cr
#' 4 PXD006617 - contains RI feature (file size 290.8 MB) \[4\]\cr\cr
#' See the Reference for further details on these examples.
# These examples are used in Bioinformatics Applications Note
# MBQN: R package for mean balanced quantile-normalization, Bioinformatics, 2018.
# @return A matrix of median- or mean-balanced quantile normalized data
#' @concept quantile, quantile normalization, rank invariance
# #' @family mbqn
#' @references \[1\] Ramond E, et al. Importance of host cell arginine uptake in Francisella phagosomal escape and ribosomal protein amounts. Mol Cell Proteomics. 2015 14(4):870-881\cr
#' \[2\] Beaumont V, et al. Phosphodiesterase 10A Inhibition Improves Cortico-Basal Ganglia Function in Huntington's Disease Models. Neuron. 2016 92(6):1220-1237\cr
#' \[3\] Turetschek R, Desalegn G, Epple T, Kaul HP, Wienkoop S. Key metabolic traits of Pisum sativum maintain cell vitality during Didymella pinodes infection: cultivar resistance and the microsymbionts' influence. J Proteomics. 2017 Mar 4. pii: S1874-3919(17)30075-1\cr
#' \[4\] Ranjbar Sistani N, Kaul HP, Desalegn G, Wienkoop S. Rhizobium Impacts on Seed Productivity, Quality, and Protection of Pisum sativum upon Disease Stress Caused by Didymella pinodes: Phenotypic, Proteomic, and Metabolomic Traits. Front Plant Sci. 2017 8:1961\cr
#' \[5\] Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note, 2018\cr
#' @examples ## Check LFQ intensities of proteomics data PXD001584 downloaded
#' ## from PRIDE for RI and NRI features
#'\dontrun{
#' mbqn.example1(which.example = 3)
#'}
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' 2017
#' @export mbqn.example1
mbqn.example1 <- function(which.example = 1, source.path = NULL){
  if(is.null(which.example)) stop("Error: Select an example between 1-4")

  ###################################################################################
  # Download proteinGroups.txt file if not already present

  #if(is.null(source.path)) {source.path = file.path(find.package("MBQN"), "examples")}
  if(is.null(source.path)) {source.path = file.path(getwd(),"examples")}

  # PXD repositories
  ids <- c("PXD001584","PXD005138","PXD005861","PXD006617")

  file <- file.path(source.path,ids[which.example],"proteinGroups.txt")
  if(!file.exists(file)){
    print("File does not exist - start download of file. This can take a few minutes...!")
    px <- rpx::PXDataset(ids[which.example])
    # General Informations on PXD File
    rpx::pxtax(px) #
    rpx::pxurl(px) # url
    rpx::pxref(px) # reference

    # list files in repository
    f <- rpx::pxfiles(px)
    ind <- grep('MaxQuant',f)

    ifelse(!dir.exists(file.path(source.path,ids[which.example])), dir.create(file.path(source.path,ids[which.example]),recursive = TRUE),'')

    if(!file.exists(f[ind])){
      fnm <- rpx::pxget(px, f[ind])
    }
    if(file.exists(f[ind])){ # unzip file
      if(length(grep("tar.gz",f[ind])>0)){
        untar(fnm, files = untar(fnm, list = T)[grepl("proteinGroups.txt",untar(fnm, list = T))], exdir=file.path(source.path,ids[which.example]))
        filesstrings::file.move(file.path(source.path,ids[which.example],untar(fnm, list = T)[grepl("proteinGroups.txt",untar(fnm, list = T))]), file.path(source.path,ids[which.example]))
        unlink(list.dirs(file.path(source.path,ids[which.example]), recursive = F), recursive = T) # delete extra folder
        unlink(fnm) # delete large .zip file
      }else{
      unzip(fnm,files = "proteinGroups.txt", exdir=file.path(source.path,ids[which.example]))
      unlink(fnm) # delete large .zip file
      }
    }
  }

  str <- unlist(strsplit(file,"/"))
  pxdid <- str[pmatch("PXD",str)]

  # Read file
  dat <- read.csv(file, allowEscapes = TRUE, check.names = FALSE,sep = "\t")

  ###################################################################################
  # Select all columns with LFQ intensities
  mtx <- as.matrix(dat[, grepl("^LFQ", names(dat))])
  mtx[mtx == 0] <- NA
  ix <- c(1:dim(mtx)[2])
  # select columns
  if(pxdid == "PXD005861") {
    ix <- c(1:9, 19:27)
   # ix <- c(10:18,31:39)
   # ix <- c(1:9, 10:18, 19:27, 31:39)
    mtx <- mtx[,ix]
    colstr <- sapply(strsplit(split = "LFQ intensity ",colnames(mtx)),function(x) x[2])
    colstr <- sub(pattern = "_20_2h",replacement = "",colstr)
    colnames(mtx) <- colstr
  }

  if(which.example==1) ylim <- c(21,36)

  # check for all proteins with NA intensity across all samples
  allColNA <- as.vector(apply(mtx, 1, function(r) {
    return(sum(!is.na(r)) == 0)
  }))
  print(paste("Number of rows with empty entries:", length(which(allColNA))))

  # check for potential contaminant and only identfied by side proteins
  ixs <- which(dat$`Potential contaminant`=="+")
  ixs <- c(ixs,which(dat$`Only identified by site`=="+"))
  ixs <- c(ixs, which(dat$Reverse=="+"))

  row.names(mtx) <- dat[, "Protein IDs"]
  featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"],
                                   proteinDescription = dat[, "Fasta headers"],
                                   #idScore = dat[, "Q-value"],
                                   isDecoy = dat[, "Reverse"] == "+",
                                   nbPeptides = dat[, "Peptides"],
                                   isNormAnchor = rep(TRUE, nrow(mtx)),
                                   isFiltered = dat[, "Reverse"] == "+",
                                   row.names = dat[, "Protein IDs"])
  ##isPotential.contaminant = dat[,"Potential contaminant"]=="+",
  featureAnnotations <- featureAnnotations[!allColNA, ]
  #mtx <- as.matrix(mtx[!allColNA, rownames(expDesign)])

  # remove empty rows
  mtx <- as.matrix(mtx[!allColNA, ]) # mod. by A. Schad

  # log2 transform intensities
  mtx <- log2(mtx)

  res <- MBQN::mbqn.check_saturation(mtx,
                                  FUN = median,
                                  show_fig = TRUE,
                                  low_thr = 0.5,
                                  filename = "",
                                  feature_index = NULL, save_fig = FALSE)

  mb <- MBQN::mbqn(as.matrix(mtx),FUN = median, na.rm = TRUE)


  # get protein name of strongest nri/ri feature
  nri_max <-as.numeric(names(which.max(res$nri)))

  featureAnnotations$proteinDescription[nri_max]
  featureAnnotations$proteinName[nri_max]
  featureAnnotations$nbPeptides[nri_max]

  df <- lapply(c(1:dim(featureAnnotations)[2]),function(j) featureAnnotations[[j]][[nri_max]])
  names(df)<- names(featureAnnotations)

  # truncate long name strings
  df$proteinName <-  paste(strtrim(df$proteinName,80),"...")
  df$proteinDescription <- paste(strtrim(df$proteinDescription,80),"...")

  # use mbqn.boxplot function, highlight RI/NRI Features
  # Create supplement Fig. S2
  save.fig <- T

  fig1.name <- paste0("fig_qnLFQ_", pxdid,".pdf")
  fig2.name <- paste0("fig_mbqnLFQ_", pxdid,".pdf")

  mbqn.mtx <- mbqn(mtx,FUN = median)
  qn.mtx <- mbqn(mtx,FUN = NULL)
  if(length(ylim)==0) ylim <- c(floor(min(range(mbqn.mtx, na.rm = T))),ceiling(max(range(mbqn.mtx, na.rm = T))))
  colnames(qn.mtx) <- colnames(mbqn.mtx) <- ix

  plot.new()
  frame()
  par(mfrow = c(1,1))
  mbqn.boxplot(mtx = qn.mtx,main = paste("QN"), ylab = "LFQ intensity", ylim = ylim,irow = c(as.numeric(names(res$nri))))
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig1.name),width=8,height=8,out.type = "pdf")
    print(paste("Save figure to ",fig1.name))
  }

# plot.new()
#  frame()
#  par(mfrow = c(1,1))
#  mbqn.boxplot(mtx = mbqn.mtx,main = paste("MBQN"), ylab = "", ylim = ylim)
#  if(save.fig){
#    dev.copy2pdf(file=file.path(getwd(),fig2.name),width=8,height=7,out.type = "pdf")
#    print(paste("Save figure to ",fig2.name))
#  }

  if(which.example==1){
  fig3.name <- paste0("fig_mbqn_vs_qn_", pxdid,".pdf")
  plot.new()
  frame()
  par(mfrow = c(1,1))
  mbqn.boxplot(mtx = mbqn.mtx,main = "MBQN with QN-feature", irow = c(as.numeric(names(res$nri)),144), vals = data.frame(QN.feature.114 = qn.mtx[144,]),ylim = ylim)
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig3.name),width=8,height=7,out.type = "pdf")
    print(paste("Save figure to ",fig3.name))
  }
  }

  # QN of data and balance only NRI/RI features
  fig4.name <- paste0("fig_qn_balance_nri_only_", pxdid,".pdf")
  plot.new()
  frame()
  par(mfrow = c(1,1))
  mbqn.boxplot(mbqn.nri(mtx,FUN = median,low_thr = 0.4),as.numeric(names(res$nri)),
               ylim = ylim, xlab= "", las=2, main = "QN with RI/NRI balanced")
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig4.name),width=8,height=7,out.type = "pdf")
    print(paste("Save figure to ",fig4.name))
  }
}
