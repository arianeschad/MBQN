#' Load and preprocess LFQ intensities
#'
#' @description Auxilary function used to load and preprocess LFQ
#' intensities of selected dataset with a PRIDE \[1\] identifier.
#' @inheritParams get_pxdfile
#' @importFrom utils read.csv
#' @details Load proteinGroup.txt file, select all samples with LFQ intensities, remove
#' empty protein features. Apply log2 transform to intensities. This function acquires source code
#' of SafeQuant::parseMaxQuantProteinGroupTxt \[2\].
#' @return List with
#' \item{\code{mtx}}{data matrix}
#' \item{\code{pxdid}}{PRIDE identifier}
#' \item{\code{featureAnnotations}}{dataframe collecting feature annotations, e.g. protein name,
#' "Potential contaminant", etc..}
#' @family example
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I, Mayer G,
#' Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.\cr
#' \[2\] See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.1.,
#' https://CRAN.R-project.org/package=SafeQuant.\cr
#' @examples ## Load LFQ intensities of proteomics data PXD001584 downloaded
#' available at PRIDE:
#'\dontrun{
#' loadFile(pxd_id = "PXD001584")
#'}
#' @author Ariane Schad
# 2017

loadFile <- function(pxd_id, source.path = NULL, file.pattern = "proteinGroups"){

  if(is.null(source.path)) {source.path = file.path(getwd())}

  file <- file.path(source.path,pxd_id,file.pattern)

  if(!file.exists(file)){
    r.input <- readline(paste("File does not exist - start downloading file. This can take a few minutes!","Do you want to continue? [y/n]",sep = "\n"))
    stopifnot(r.input=="y")
    print("Proceed with download...")
    ifelse(!dir.exists(file.path(source.path,pxd_id)), dir.create(file.path(source.path,pxd_id),recursive = TRUE),'')
    # Download proteinGroups.txt file if not already present
    get_pxdfile(pxd_id = pxd_id, source.path = source.path, file.pattern = file.pattern)

  } # fi file.exist
  ##################################################################################
  str <- unlist(strsplit(file,"/"))
  pxdid <- str[pmatch("PXD",str)]

  # Read file
  dat <- read.csv(file, allowEscapes = TRUE, check.names = FALSE,sep = "\t")

  # Select all columns with LFQ intensities
  mtx <- as.matrix(dat[, grepl("^LFQ", names(dat))])
  mtx[mtx == 0] <- NA
  ix <- c(1:dim(mtx)[2])
  ylim <- NULL

  # check for proteins with NA intensity across all samples
  allColNA <- as.vector(apply(mtx, 1, function(r) {
    return(sum(!is.na(r)) == 0)
  }))
  print(paste("Number of proteins with empty entries:", length(which(allColNA))))

  # check for potential contaminant and only identfied by side proteins
  ixs <- which(dat$`Potential contaminant`=="+")
  ixs <- c(ixs,which(dat$`Only identified by site`=="+"))
  ixs <- c(ixs, which(dat$Reverse=="+"))

  featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"],
                                   proteinDescription = dat[, "Fasta headers"],
                                   #idScore = dat[, "Q-value"],
                                   isDecoy = dat[, "Reverse"] == "+",
                                   nbPeptides = dat[, "Peptides"],
                                   isNormAnchor = rep(TRUE, nrow(mtx)),
                                   isFiltered = dat[, "Reverse"] == "+",
                                   isPotential.contaminant = dat[,"Potential contaminant"]=="+")
                                   #row.names = dat[, "Protein IDs"])

  featureAnnotations <- featureAnnotations[!allColNA, ]

  row.names(mtx) <- dat[, "Protein IDs"]

  # remove empty rows
  mtx <- as.matrix(mtx[!allColNA, ])

  # log2 transform intensities
  mtx <- log2(mtx)

  return(list(mtx = mtx, pxdid = pxdid, featureAnnotations = featureAnnotations))
}

