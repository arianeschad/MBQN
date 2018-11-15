#' Example illustrating the analysis of PRIDE data with the MBQN package
#'
#' @description This function runs an example to illustrate the check of LFQ intensities
#' of proteins for rank invariant (RI) and nearly rank invariant (NRI) features
#' and the normalization of the data. It downloads a selected dataset from PRIDE
#' or it uses the dataset in filepath. The example needs the R package "rpx" to
#' download data from PRIDE.
#' @param which.example Numerical values between 1-4, default = 1,
#' to select between 4 dataset.
#' @param filepath Character with path where to store and search for proteinGroups.txt
#' file; default = NULL uses "installationpath/MBQN/examples".
#' @inheritParams mbqn
#' @importFrom grDevices dev.copy2pdf dev.off dev.size pdf
#' @details Collecting information on the experiments requires the package rpx by Laurent Gatto (2017),
# #' rpx: R Interface to the ProteomeXchange Repository.
# #' R package version
#' version 1.10.2 from https://github.com/lgatto/rpx.\cr
#' The function used to read the proteinGroups.txt files uses source code
#' of SafeQuant::parseMaxQuantProteinGroupTxt, \cr
# #' See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.1.
#' version 2.3.1, https://CRAN.R-project.org/package=SafeQuant.\cr\cr
#' Experiments on PRIDE that are supported by this Example:\cr
#' 1 PXD001584 - contains a RI (default)\cr
#' 2 PXD005138 - contains a RI feature\cr
#' 3 PXD005861 - contains RI feature\cr
#' 4 PXD006617 - contains RI feature\cr\cr
#' See the Reference for further details on these examples.
# These examples are used in Bioinformatics Applications Note
# MBQN: R package for mean balanced quantile-normalization, Bioinformatics, 2018.
# @return A matrix of median- or mean-balanced quantile normalized data
#' @concept quantile, quantile normalization, rank invariance
#' @family mbqn
# #' @include mbqn
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. Bioinf. Appl. Note, 2018
#'@examples ## Check LFQ intensities of proteomics data PXD001584 downloaded
#'from PRIDE for RI and NRI features
#'\dontrun{
#' example1(which.example = 1)
#'}
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' 2017
#' @export mbqn.example1
# Example: Analyse selected PRIDE data used in the Bioinformatics Applications Note
# "MBQN: R package for mean balanced quantile-normalization, Bioinformatics, 2018"
# This script depends on R packege "rpx" and "SafeQuant"
#
# Collecting information on the experiments requires the rpx package
# install.packages("rpx")
# Laurent Gatto (2017). rpx: R Interface to the ProteomeXchange Repository. R package version
# 1.10.2. https://github.com/lgatto/rpx
#
# The function used to read the proteinGroups.txt files is based on parts of the source code
# of SafeQuant::parseMaxQuantProteinGroupTxt
# See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.1.
# https://CRAN.R-project.org/package=SafeQuant
# Usage:
# > source("examples/example1.R")
# > example1(which.example = 1)
# Arguments:
# which.example - select an example between 1-4, default = 1
# filepath - path where file is stored and found, default = NULL specifies "installationpath/MBQN/examples"
# PXD experiments selected as examples
# 1. PXD001584 - contains a RI
# file = '/Users/schad/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2015/01/PXD001584/MaxQuantOutput/proteinGroups.txt'
# 2. PXD005138 - contains a RI feature
# 3. PXD005861 - contains RI feature (default)
# 4. PXD006617 - contains RI feature
#file = '~/examples/PXD001584/proteinGroups.txt'
#file = '~/examples/PXD005138/proteinGroups.txt'
#file = '~/examples/PXD005861/proteinGroups.txt'
#file = '~/examples/PXD006617/proteinGroups.txt'
# A. Schad, April 2018
mbqn.example1 <- function(which.example = 1, source.path = NULL){
  if(is.null(which.example)) stop("Error: Select an example between 1-4")

  ###################################################################################
  # Download proteinGroups.txt file if not already present

  if(is.null(source.path)) {source.path = file.path(find.package("MBQN"), "examples")}

  # PXD repositories
  ids <- c("PXD001584","PXD005138","PXD005861","PXD006617")

  file <- file.path(source.path,ids[which.example],"proteinGroups.txt")
  if(!file.exists(file)){

    px <- rpx::PXDataset(ids[which.example])
    # General Informations on PXD File
    rpx::pxtax(px)
    rpx::pxurl(px)
    rpx::pxref(px)

    # list files in repository
    f <- rpx::pxfiles(px)
    ind <- grep('MaxQuant',f)

    ifelse(!dir.exists(file.path(source.path,"examples",ids[which.example])), dir.create(file.path(source.path,"examples",ids[which.example])),'')

    if(!file.exists(f[ind])){
      fnm <- rpx::pxget(px, f[ind])
    }
    if(file.exists(f[ind])){ # unzip file
      unzip(fnm,files = "proteinGroups.txt", exdir=file.path(source.path,"examples",ids[which.example]))
      unlink(fnm) # delete large .zip file
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

  # select columns
  if(pxdid == "PXD001584") {
    ix <- c(1:9, 19:27)
  }

  res <- MBQN::mbqn.check_saturation(mtx[,ix],
                                  FUN = median,
                                  show_fig = TRUE,
                                  low_thr = 0.5,
                                  filename = "",
                                  feature_index = NULL, save.fig = FALSE)

  mb <- MBQN::mbqn(as.matrix(mtx[,ix]),FUN = median, na.rm = TRUE)


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

  # use mbqn.boxplot function, highlight RI Feature

  # Create supplement Fig. S2
  save.fig <- T

  fig1.name <- paste0("fig_qnLFQ_", pxdid,".pdf")
  fig2.name <- paste0("fig_mbqnLFQ_", pxdid,".pdf")

  mbqn.mtx <- mbqn(mtx[,ix],FUN = median)
  qn.mtx <- mbqn(mtx[,ix],FUN = NULL)
  colnames(qn.mtx) <- colnames(mbqn.mtx) <- ix

  plot.new()
  frame()
  par(mfrow = c(1,1))
  mbqn.boxplot(mtx = qn.mtx,main = paste("QN"), ylab = "LFQ intensity")
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig1.name),width=8,height=7,out.type = "pdf")
    print(paste("Save figure to ",fig1.name))
  }
  #mbqn.boxplot(mtx = mtx[,ix],main = paste("MBQN of LFQ intensities for", pxdid))
  plot.new()
  frame()
  par(mfrow = c(1,1))
  mbqn.boxplot(mtx = mbqn.mtx,main = paste("MBQN"), ylab = "")
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig2.name),width=8,height=7,out.type = "pdf")
    print(paste("Save figure to ",fig2.name))
  }


  fig3.name <- paste0("fig_mbqn_vs_qn_", pxdid,".pdf")
  plot.new()
  frame()
  par(mfrow = c(1,1))
  mbqn.boxplot(mtx = mbqn.mtx,main = paste("MBQN vs. QN"), irow = c(res$var0_feature,144,600), vals = qn.mtx[144,])
  if(save.fig){
    dev.copy2pdf(file=file.path(getwd(),fig3.name),width=8,height=7,out.type = "pdf")
    print(paste("Save figure to ",fig3.name))
  }

  # standard deviation
  breaks = seq(0,2,0.01)
  hist(apply(mbqn.mtx,1,sd, na.rm = TRUE),main = "Histogram of feature variation",
       xlab = "std",breaks= breaks,col = 3,add = F)
  hist(apply(qn.mtx,1,sd, na.rm = TRUE),breaks = breaks,col = 4,add = TRUE)
  legend(x = 1.5, y = 50, legend = c("MBQN","QN"),fill = c(3,4), bty = "n",border = FALSE)

}
