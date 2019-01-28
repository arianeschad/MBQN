## Subfunction to load selected PRIDE dataset and preprocess it, i.e.,
# load proteinGroup intensity data, select LFQ intensities, remove
# empty features, and log2 transform intensities

mbqnLoadExample <- function(which.example = 1, source.path = NULL){

  if(is.null(which.example)) stop("Error: Select an example between 1-5")

  # PXD repositories
  ids <- c("PXD001584","PXD005138","PXD005861","PXD006617")

  ########################################################################
  # Download proteinGroups.txt file if not already present

  if(is.null(source.path)) {source.path = file.path(getwd())}

  file <- file.path(source.path,ids[which.example],"proteinGroups.txt")
  if(!file.exists(file)){
    r.input <- readline(paste("File does not exist - start downloading file. This can take a few minutes!","Do you want to continue? [y/n]",sep = "\n"))
    stopifnot(r.input=="y")
    print("Proceed with download...")

    px <- rpx::PXDataset(ids[which.example])

    # General Informations on PXD File
    rpx::pxtax(px)
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
  } # fi file.exist

  str <- unlist(strsplit(file,"/"))
  pxdid <- str[pmatch("PXD",str)]

  # Read file
  dat <- read.csv(file, allowEscapes = TRUE, check.names = FALSE,sep = "\t")

  ###################################################################################
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
                                   row.names = dat[, "Protein IDs"])
  # isPotential.contaminant = dat[,"Potential contaminant"]=="+",
  featureAnnotations <- featureAnnotations[!allColNA, ]

  row.names(mtx) <- dat[, "Protein IDs"]

  # remove empty rows
  mtx <- as.matrix(mtx[!allColNA, ])

  # log2 transform intensities
  mtx <- log2(mtx)

  return(list(mtx = mtx, featureAnnotations = featureAnnotations))
}

