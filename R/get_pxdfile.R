#' Download data from PRIDE
#'
#' @description Auxiliary function which downloads selected data, like protein group intensities, from PRIDE \[1\].
#' @param pxd_id character specifying the PRIDE identifier.
#' @param source.path character specifying the path where to store and search for the
#' "proteinGroups.txt"-file; default = NULL uses current working directory.
#' @param file.pattern character specifying the kind of dataset for download, e.g.
#' "proteinGroups.txt" or "peptides.txt".
#' @importFrom filesstrings file.move
#' @importFrom utils untar unzip
#' @details This function requires the R package rpx \[2\].
# #' The function used to read the proteinGroups.txt files acquires source code
# #' of SafeQuant::parseMaxQuantProteinGroupTxt, \cr
# #' See "SafeQuant" by Erik Ahrne (2016). SafeQuant: A Toolbox for the Analysis of Proteomics Data. R package version 2.3.1.
# #' version 2.3.1, https://CRAN.R-project.org/package=SafeQuant.\cr\cr
#' @references
#' \[1\] Vizcaíno JA, Csordas A, del-Toro N, Dianes JA, Griss J, Lavidas I, Mayer G,
#' Perez-Riverol Y, Reisinger F, Ternent T, Xu QW, Wang R, Hermjakob H.
#' 2016 update of the PRIDE database and related tools. Nucleic Acids Res.
#' 2016 Jan 1;44(D1): D447-D456. PubMed PMID:26527722.\cr
#' \[2\] See "rpx" by Laurent Gatto (2017). rpx: R Interface to the ProteomeXchange Repository.
#' Package version 1.10.2, https://github.com/lgatto/rpx.
#' @examples ## Download protein LFQ intensities of PXD001584 from PRIDE
#'\dontrun{
#' get_pxdfile(pxd_id = "PXD001584")
#'}
#' @author Ariane Schad
# 2018
get_pxdfile <- function(pxd_id, source.path = NULL, file.pattern = "proteinGroups"){
  px <- rpx::PXDataset(pxd_id)

  if(is.null(source.path)) {source.path = file.path(getwd())}

  # Check if package preprocessCore is installed  to run this function
  if (!requireNamespace("rpx", quietly = TRUE)) {
      stop("Package \"pkg\" is required for this function to work. Please install it.",
           call. = FALSE)}

  # General Informations on PXD File
  rpx::pxtax(px)
  rpx::pxurl(px) # url
  rpx::pxref(px) # reference

  # list files in repository
  f <- rpx::pxfiles(px)
  ind <- grep('MaxQuant',f)

  if(!file.exists(f[ind])){
    fnm <- rpx::pxget(px, f[ind])
  }
  if(file.exists(f[ind])){ # unzip file
    if(length(grep("tar.gz",f[ind])>0)){
      files = untar(fnm, list = TRUE)[grepl(file.pattern,untar(fnm, list = TRUE))]
      untar(fnm,
            files = files,
            exdir=file.path(source.path,pxd_id))
      filesstrings::file.move(file.path(source.path,pxd_id, files),
                              file.path(source.path,pxd_id))
      unlink(list.dirs(file.path(source.path,pxd_id), recursive = FALSE), recursive = TRUE) # delete extra folder
      unlink(fnm) # delete large .zip file
    }else{ #.gz
     # files = unzip(fnm, list = TRUE)[grepl(file.pattern,unzip(fnm, list = TRUE))]
      unzip(fnm,files = "proteinGroups.txt", exdir=file.path(source.path,pxd_id))
      unlink(fnm) # delete large .zip file
    }
  }

}
