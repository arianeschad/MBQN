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
#' @return Datafile from PRIDE.
#' @details This function requires the R package rpx \[2\].
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I,
#' Mayer G, Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.\cr
#' \[2\] See "rpx" by Laurent Gatto (2017). rpx: R Interface to the
#' ProteomeXchange Repository. Package version 1.10.2,
#' https://github.com/lgatto/rpx.
#' @examples ## Download protein LFQ intensities of PXD001584 from PRIDE
#'\dontrun{
#' getPXDfile(pxd_id = "PXD001584")
#'}
#' @author Ariane Schad
# 2018
getPXDfile <- function(pxd_id, source.path = NULL,
                        file.pattern = "proteinGroups"){
  px <- rpx::PXDataset(pxd_id)

  if(is.null(source.path)) {source.path = file.path(getwd())}

  # Check if package preprocessCore is installed  to run this function
  if (!requireNamespace("rpx", quietly = TRUE)) {
      stop("Package \"pkg\" is required for this function to work.
           Please install it this package first!", call. = FALSE)}

  # General Informations on PXD File
  rpx::pxtax(px)
  rpx::pxurl(px) # url
  rpx::pxref(px) # reference

  # list files in repository
  f <- rpx::pxfiles(px)
  ind <- grep('MaxQuant',f)

  if(!file.exists(f[ind])){ # if it does not yet exist locally
    fnm <- rpx::pxget(px, f[ind], method="libcurl") # download
  }
  if(file.exists(f[ind])){ # if zip file is available, start unzipping
    if(length(grep("tar.gz",f[ind])>0)){
      files = untar(fnm, list = TRUE)[grepl(file.pattern,untar(fnm,
                                                 list = TRUE))]
      filepath = file.path(source.path,pxd_id)
      untar(fnm,
            files = files,
            exdir=)
      file.rename(file.path(source.path,files),
                  file.path(filepath,"proteinGroups.txt"))
      unlink(list.dirs(filepath, recursive = FALSE),
             recursive = TRUE) # delete extra folder
      unlink(fnm) # delete large .zip file
    }else{ #.gz
      unzip(fnm,files = "proteinGroups.txt",
            exdir=filepath)
      unlink(fnm) # delete large .zip file
    }
  }else{
    message("No MaxQuant file available.")
  }

}
