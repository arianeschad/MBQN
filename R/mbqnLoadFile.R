#' Load and preprocess LFQ intensities
#'
#' @description Auxilary function used to load and preprocess LFQ
#' intensities of selected dataset with a PRIDE \[1\] identifier.
#' @inheritParams getPXDfile
#' @importFrom utils read.csv
#' @details Load proteinGroup.txt file, select all samples with LFQ intensities, remove
#' empty protein features. Apply log2 transform to intensities. This function acquires source code
#' of SafeQuant::parseMaxQuantProteinGroupTxt \[2\].
#' @return List with
#' \item{\code{mtx}}{data matrix}
#' \item{\code{pxdid}}{PRIDE identifier}
#' \item{\code{featureAnnotations}}{dataframe collecting feature annotations, e.g. protein name,
#' "Potential contaminant", etc..}
#' \item{\code{ixs}}{locical array indicating rows of potential contaminant features, or features identified only by site, reverse.}
#' @seealso [getPXDfile()] for downloading data from PRIDE.
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I, Mayer G,
#' Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.\cr
#' \[2\] See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the Analysis
#' of Proteomics Data. R package version 2.3.1.,
#' https://CRAN.R-project.org/package=SafeQuant.\cr
#' @examples ## Load LFQ intensities of proteomics data of PXD001584:
#' \dontrun{
#' mbqnLoadFile(pxd_id = "PXD001584")
#' }
#' @author Ariane Schad
#  2017
#' @export mbqnLoadFile
mbqnLoadFile <- function(pxd_id, source.path = NULL, file.pattern = "proteinGroups.txt"){

  if(is.null(source.path)) {source.path = file.path(getwd())}

  fdir <- file.path(source.path,pxd_id)
  file <- list.files(fdir, pattern = file.pattern,full.names = TRUE)

  if(length(file)==0){ # file does not exist in fdir
    r.input <- readline(paste("File does not exist - start downloading file. This can take a few minutes!","\n",
                              "Do you want to continue? [y/n]", "\n"))
    stopifnot(r.input=="y")
    message("Proceed with download...")
    ifelse(!dir.exists(fdir), dir.create(fdir,recursive = TRUE),'')
    # Download proteinGroups.txt file if not already present
    getPXDfile(pxd_id = pxd_id, source.path = source.path, file.pattern = file.pattern)
    file <- list.files(fdir, pattern = file.pattern,full.names = TRUE)
  } # fi file.exist
  ##################################################################################
  str <- unlist(strsplit(file,"/"))
  pxdid <- str[pmatch("PXD",str)]

  # Read file
  dat <- read.csv(file, allowEscapes = TRUE, check.names = FALSE,sep = "\t")

  # Select all columns with LFQ intensities
  mtx <- as.matrix(dat[, grepl("^LFQ", names(dat))])
  mtx[mtx == 0] <- NA

  # check for proteins with NA intensity across all samples
  allColNA <- as.vector(apply(mtx, 1, function(r) {
    return(all(is.na(r)))
  }))
  message(paste("Number of proteins with empty entries:", length(which(allColNA))))


  # check if exist and append to featureAnnotations
  featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"])
  annotations <- c("Fasta headers", "Q-value", "Reverse", "Peptides", "Reverse",
                   paste(c("Potential contaminant","Contaminant"),collapse = "|"),
                   "Only identified by site")
  fieldnames <- c("proteinDescription", "idScore","isDecoy", "nbPeptides",
                  "isFiltered", "isPotential.contaminant", "isIdentified.by.site")


  # check for potential contaminant and only identfied by side proteins
  bool1 <- dat[,grep(annotations[6],colnames(dat),value = TRUE)]=="+"
  bool2 <- dat[["Only identified by site"]]=="+" & !is.na(dat[["Only identified by site"]])
  bool3 <- dat[["Reverse"]]=="+"
  ixs <-  bool1 | bool2 | bool3  # logical is better because it has same dimension as the data

  for (i in seq_len(length(fieldnames))){
    if(length(grep(annotations[i],colnames(dat)))>0){
      featureAnnotations[[fieldnames[i]]] <- ifelse(length(grep(fieldnames[i],
                                                                c("isDecoy","isPotential.contaminant",
                                                                  "Only identified by site")))>0,
                                                   dat[,grep(annotations[i],colnames(dat),value = TRUE)]=="+",
                                                   dat[,annotations[i]])
       }
  }


  featureAnnotations <- featureAnnotations[!allColNA, ]

  row.names(mtx) <- dat[, "Protein IDs"]

   # remove empty rows
  mtx <- as.matrix(mtx[!allColNA, ])

  # log2 transform intensities
  mtx <- log2(mtx)

  
  # SummarizedExperiment(assays=list(intensities=mtx),
  #                      row.names=featureAnnotations, colData=colData)
  # sExp <- SummarizedExperiment(mtx, featureNames=featureAnnotations, sampleNames=ixs)
  # adata
  return(list(mtx = mtx, pxdid = pxdid, featureAnnotations = featureAnnotations, ixs = ixs))
}

