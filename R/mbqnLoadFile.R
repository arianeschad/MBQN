#' Load and preprocess LFQ intensities
#'
#' @description Auxilary function used to load and preprocess LFQ
#' intensities of selected dataset with a PRIDE \[1\] identifier.
#' @inheritParams getPXDfile
#' @importFrom utils read.csv
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @details Load proteinGroup.txt file, select all samples with LFQ
#' intensities, remove empty protein features. Apply log2 transform to
#' intensities. This function acquires source code of
#' SafeQuant::parseMaxQuantProteinGroupTxt \[2\].
#' @return SummarizedExperiment with
#' \item{\code{data}}{data matrix}
#' \item{\code{pxdid}}{PRIDE identifier}
#' \item{\code{featureAnnotations}}{dataframe collecting feature annotations,
#' e.g. protein name, "Potential contaminant", etc..}
#' \item{\code{ixs}}{locical array indicating rows of potential contaminant
#' features, or features identified only by site, reverse.}
#' @seealso [getPXDfile()] for downloading data from PRIDE.
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I,
#' Mayer G, Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.\cr
#' \[2\] See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the
#' Analysis of Proteomics Data. R package version 2.3.1.,
#' https://CRAN.R-project.org/package=SafeQuant.\cr
#' @examples ## Load LFQ intensities of proteomics data of PXD001584:
#' library(SummarizedExperiment)
#' out <- mbqnLoadFile(pxd_id = "PXD001584")
#' mbqn(assays(out)[["data"]])
#' @author Ariane Schad
#  2017
#' @export mbqnLoadFile
mbqnLoadFile <- function(
    pxd_id, source.path = NULL, file.pattern = "proteingroups"){
    
    if (is.null(source.path)) {source.path = file.path(getwd())}
    
    fdir <- file.path(source.path,pxd_id)
    file <- list.files(
        fdir, pattern = file.pattern,full.names = TRUE, recursive= TRUE, ignore.case = TRUE)
    
    if (length(file)==0){ # file does not exist in fdir
        message("File does not exist - proceed with download...")

        # Download proteinGroups.txt file if not already present
        status <- getPXDfile(pxd_id = pxd_id, source.path = source.path,
                             file.pattern = file.pattern)
        file <- list.files(
            fdir, pattern = file.pattern,full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
    } # fi file.exist
    ############################################################################
    str <- unlist(strsplit(file,"\""))
    if (length(str)==1){
        str <- unlist(strsplit(file,"/"))
    }
    pxdid <- str[pmatch("PXD",str)]
    
    # Read file
    dat <- read.csv(file[1], allowEscapes = TRUE, check.names = FALSE,sep = "\t")
    
    # Select all columns with LFQ intensities
    mtx <- as.matrix(dat[, grepl("^LFQ", names(dat))])
    mtx[mtx == 0] <- NA
    mtx[mtx == "NaN"] <- NA
    # convert character matrix to numeric matrix, strings that are possibly 
    # left in matrix by accident are converted to NAs
    mtx<- `dimnames<-`(`dim<-`(as.numeric(mtx), dim(mtx)), dimnames(mtx))
    
    # check for proteins with NA intensity across all samples
    allColNA <- as.vector(apply(mtx, 1, function(r) {
        return(all(is.na(r)))
    }))
    message(paste("Number of proteins with empty entries:",
                  length(which(allColNA))))
    
    if (!("Protein IDs" %in% colnames(dat))){
        colnames(dat)[1] <- "Protein IDs"
    }
    # check if exist and append to featureAnnotations
    featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"])
    annotations <- c(
        "Fasta headers", "Q-value", "Reverse", "Peptides", "Reverse",
        paste(c("Potential contaminant","Contaminant"),collapse = "|"),
        "Only identified by site")
    fieldnames <- c("proteinDescription", "idScore","isDecoy", "nbPeptides",
                    "isFiltered", "isPotential.contaminant", "isIdentified.by.site")
    
    # check for potential contaminant and only identfied by side proteins
    bool1 <- dat[,grep(annotations[6],colnames(dat),value = TRUE)]=="+" & !is.na(
        dat[,grep(annotations[6],colnames(dat),value = TRUE)])
    bool2 <- dat[["Only identified by site"]]=="+" & !is.na(
        dat[["Only identified by site"]])
    bool3 <- dat[["Reverse"]]=="+" & !is.na(dat[["Reverse"]])
    # logical is better because it has same dimension as the data
    ixs <-  bool1 | bool2 | bool3
    
    for (i in seq_len(length(fieldnames))){
        if(length(grep(annotations[i],colnames(dat)))>0){
            if(length(grep(fieldnames[i],c("isDecoy",
                                           "isPotential.contaminant",
                                           "Only identified by site"))) > 0){
                featureAnnotations[[fieldnames[i]]] <-
                    dat[,grep(annotations[i],colnames(dat),value = TRUE)]=="+"
            } else {
                featureAnnotations[[fieldnames[i]]] <- dat[,annotations[i]]
            }
        }
    }
    
    # featureAnnotations <- featureAnnotations[!allColNA, ]
    row.names(mtx) <- dat[, "Protein IDs"]
    
    # # remove empty rows
    # mtx <- as.matrix(mtx[!allColNA, ])
    
    # log2 transform intensities
    mtx <- log2(mtx)
    
    featureAnnotations <- as.list(featureAnnotations)
    featureAnnotations[["ixs"]]<-ixs
    
    sexp <- SummarizedExperiment(assays=list(data=mtx),
                                 rowData=featureAnnotations,
                                 metadata=list(pxdid = pxdid))
    # remove empty rows
    sexp <- sexp[!allColNA, ]
    return(sexp)
}
