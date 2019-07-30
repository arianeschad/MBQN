#' Download data from PRIDE
#'
#' @description Auxiliary function which downloads selected data, like protein
#' group intensities, from PRIDE \[1\].
#' @param pxd_id the PRIDE identifier.
#' @param source.path pathname where to store and search for the
#' data file, e.g. "proteinGroups.txt"-file; default = NULL uses current working
#' directory.
#' @param file.pattern character specifying the kind of dataset for download,
#' e.g. "proteinGroups" or "peptides".
#' @importFrom utils untar unzip
#' @importFrom BiocFileCache bfcadd
#' @return status (0=ok, 1=not MaxQuant data set, 2=no proteinGroup file)
#' @details This function requires the R packages rpx \[2\] and BiocFileChace.
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I,
#' Mayer G, Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.\cr
#' \[2\] See "rpx" by Laurent Gatto (2017). rpx: R Interface to the
#' ProteomeXchange Repository. Package version 1.10.2,
#' https://github.com/lgatto/rpx.
#' @examples ## Download protein LFQ intensities of PXD001584 from PRIDE
#' getPXDfile(pxd_id = "PXD001584")
#' @author Ariane Schad
# 2018
#' @export getPXDfile
getPXDfile <- function(pxd_id, source.path = NULL,
                       file.pattern = "proteinGroups"){
    px <- rpx::PXDataset(pxd_id)
    
    if (is.null(source.path)) {source.path = file.path(getwd())}
    pdxFolder = file.path(source.path,pxd_id)
    
    # Check if package preprocessCore is installed  to run this function
    if (!requireNamespace("rpx", quietly = TRUE)) {
        stop("Package \"pkg\" is required for this function to work.
            Please install it this package first!", call. = FALSE)}
    
    # General Informations on PXD File:
    # rpx::pxtax(px)
    # rpx::pxurl(px) # url
    # rpx::pxref(px) # reference
    
    status <- 0 # everything ok
    
    # list files in repository
    repoFiles <- rpx::pxfiles(px)
    #ind <- grep('MaxQuant',repoFiles)
    ind <- grep('^(?!.*peptide).*maxquant|proteingroups.*$', tolower(repoFiles), perl=T)
    if (length(ind)==0){
        message("Only MaxQuant data formats are supported. Download stopped.")
        status <- 1
        return(status)
    }
    if (length(ind)>1)
        ind <- ind[1] # only the first match
    
    destDirName <- "MBQN"
    destDir <- rappdirs::user_cache_dir(appname=destDirName)
    bfc <- BiocFileCache::BiocFileCache(destDir, ask=TRUE)
    destFile <- BiocFileCache::bfcadd(bfc,
                                      pxd_id,
                                      fpath=paste(rpx::pxurl(px),
                                                  repoFiles[ind],sep="/"))
    
    if (file.exists(destFile)){ # if zip file is available, start unzipping
        if(length(grep(".txt",repoFiles[ind])>0)){
            dir.create(file.path(pdxFolder), showWarnings = FALSE)
            # files = repoFiles[ind]
            file.copy(destFile, pdxFolder)
            unlink(destFile)
            
        } else if (length(grep("tar.gz",repoFiles[ind])>0)){
            files = untar(destFile, list = TRUE) # which files are in archive
            files = files[grepl(file.pattern,files)] # only pattern match
            
            untar(destFile, files = files, exdir=)
            
            # Change directory name to pxd id
            tryCatch(
                {
                    file.rename(file.path(source.path, sub("/[^/]+$", "", files)),
                                file.path(pdxFolder))
                },
                error = function(e){
                    print(e)
                }
            )
            
            
            unlink(list.dirs(pdxFolder, recursive = FALSE),
                   recursive = TRUE) # delete extra folder
            unlink(destFile) # delete large .zip file
            # delete download folder
            unlink(destDir,recursive = TRUE, force = TRUE)
        } else { #.gz
            files = unzip(destFile, list = TRUE) # which files are in archive
            # only pattern match
            files = files$Name[grepl(file.pattern,files$Name)]
            
            if (length(files)>0){
                unzip(destFile,files = files,exdir=pdxFolder)
                
                # Change directory name to pxd id
                tryCatch(
                    {
                        file.rename(file.path(source.path, sub("/[^/]+$", "", files)),
                                    file.path(pdxFolder))
                    },
                    error = function(e){
                        print(e)
                    }
                )
                
                unlink(destFile) # delete large .zip file
                # delete download folder
                unlink(destDir,recursive = TRUE, force = TRUE)
            } else {
                message(sprintf(
                    "File %s not found in the downloaded data archive.",
                    file.pattern))
                status <- 3
            }
        }
    } else {
        message("No MaxQuant file available.")
        status <- 2
    }
    return(status)
}
