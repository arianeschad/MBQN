#' Generate a random/structured data matrix
#'
#' @description Generate a random data matrix with or without proteomics-like properties
#' @param model Two different type of matrix models are available: "rand" generates a random matrix
#' of size nrow x ncol, "omics" generates a Gaussian random matrix which mimics intensity profiles and
#' missing values as present in a selected real data set.
#' @param nrow number of rows of data matrix.
#' @param ncol number of columns of data matrix.
#' @param seed Seed for random number generator. Default \code{seed <- 1234}.
#' @details For model "rand" the random matrix entries are drawn from a standard normal distribution
#' of mean 0 and sigma 1. For model "omics" the meand and standard deviation of each row i is
#' drawn from Gaussian distribution \eqn{N(\mu_i,sd_i^2)} with \eqn{sd_i~N(0,0.0625)} and \eqn{\mu_i~N(28,4)}.
#' About 21\% of the values are omitted according to the missing value pattern present in the LFQ intensities
#' of PXD001584 (Ramond et al. 2015).
#' @return \code{matrix} of size nrow x ncol
#' @importFrom utils read.csv untar unzip
#' @importFrom stats rnorm
#' @importFrom graphics image layout points rect
# #' @keywords quantile normalization proteomics
#' @references Schad, A. and Kreuz, C., MBQN: R package for mean balanced quantile normalization. In prep. 2019. \cr
#' Ramond, E. et al. (2015) Importance of host cell arginine uptake in Francisella phagosomal
#' escape and ribosomal protein amounts. Mol Cell Proteomics 14, 870-881.
#' @examples
#' mbqn.simu.data(model = "rand", 1000,10)
#' mbqn.simu.data(model = "rand")
#' mbqn.simu.data(model = "omics")
#' @author A. Schad, \email{ariane.schad@zbsa.de}
#' @export mbqn.simu.data
mbqn.simu.data <- function(model = "rand", nrow = NULL, ncol = NULL,seed = 1234){

  if(is.null(nrow)) nrow <- 1000
  if(is.null(ncol)) ncol <- 10
  if(!is.null(seed)) set.seed(seed)

  if(model=="rand"){
    # generate a random matrix without NAs
    dat <- replicate(ncol, rnorm(nrow))
  }

  if(model=="omics"){
    # generate a structured random matrix with NAs
    # the NA structure is obtained from a real dataset

    # load a real dataset from a PRIDE archive
    mtx <- mbqn.load.example()
    dat <- replicate(dim(mtx)[2], rnorm(dim(mtx)[1])*0.25)+rnorm(dim(mtx)[1])*2+28

    s <-sort(apply(dat, 1,mean, na.rm =T), index.return = T , decreasing = F)
    # image(t(dat.na[s$ix,]))
    s.mtx <-sort(apply(mtx, 1,mean, na.rm =T), index.return = T , decreasing = F)
    mtx <- mtx[s.mtx$ix,]
    dat <- dat[s$ix,]
    dat[is.na(mtx)] <- NA

    dev.off()
    plot.new()
    frame()
    par(mfrow=c(2,2))
    image(t(mtx), xlab = "sample", ylab = "feature row (sorted)", main = "real measurements")
    image(t(dat), xlab = "sample", main = "simulated")

    plot(apply(mtx,1,mean,na.rm =T),apply(is.na(mtx),1,mean,na.rm =T), xlab = "intensity", ylab = "NA probability")
    points(apply(dat,1,mean,na.rm =T),apply(is.na(dat),1,mean,na.rm =T), col =2)
    legend("bottomleft",legend = c("data","simulated"),  pch = 1, col = c(1,2), bty ="n", y.intersp = 1.5)
  }

  return(dat)

}


## subfunction to load dataset
mbqn.load.example <- function(which.example = 1, source.path = NULL){
  if(is.null(which.example)) stop("Error: Select an example between 1-4")

  ###################################################################################
  # Download proteinGroups.txt file if not already present

  if(is.null(source.path)) {source.path = file.path(getwd(),"examples")}

  # PXD repositories
  ids <- c("PXD001584","PXD005138","PXD005861","PXD006617")

  file <- file.path(source.path,ids[which.example],"proteinGroups.txt")
  if(!file.exists(file)){
    print("File does not exist - start downloading file. This can take a few minutes...!")

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
  ylim <- NULL

  # select columns
  if(pxdid == "PXD001584") {
    ix <- c(1:9, 19:27)
    mtx <- mtx[,ix]
    ylim <- c(22.5,36)
  }

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

  # remove empty rows
  mtx <- as.matrix(mtx[!allColNA, ])

  # log2 transform intensities
  mtx <- log2(mtx)
  return(mtx)
}




