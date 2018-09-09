# Examples from PRIDE data used in the Bioinformatics Applications Note
# "MBQN: R package for mean balanced quantile-normalization, Bioinformatics, 2018"

# The function to read proteinGroups.txt files is based on the source code
# of SafeQuant::parseMaxQuantProteinGroupTxt
# install.packages("SafeQuant")
# See details of function SafeQuant::parseMaxQuantProteinGroupTxt

# Collecting information on the experiments requires the rpx package
# install.packages("rpx")

# A. Schad, April 2018

rm(list = ls())
# requires rpx
# library(rpx)
# ids <- c("PXD001584","PXD005138","PXD005861","PXD006617")
#
 # for (i in 1:length(ids)){
 #   px <- rpx::PXDataset(ids[i])
 #   rpx::pxtax(px)
 #   rpx::pxurl(px)
 #   rpx::pxref(px)
 # }


####################
# pxfiles(px)
#fnm <- pxget(px, "README.txt")
#fnm <- pxget(px, "PXD000001_mztab.txt")

####################

#library("MSnbase")
#readMzTabData(fnm, "PEP")
#unlink("PXD000001_mztab.txt")

# read.table(file.path(""), header = TRUE)
##################################

#install.packages("SafeQuant")
#load("SafeQuant")

#see details of function
#SafeQuant::parseMaxQuantProteinGroupTxt

#############

# contains a RI
#file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2015/01/PXD001584/MaxQuantOutput/proteinGroups.txt'

# contains a RI feature, interesting file
#file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2016/12/PXD005138/proteinGroups.txt'

# contains RI feature, interessantes file!
file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2017/03/PXD005861/proteinGroups.txt'

# contains RI feature
#file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2017/12/PXD006617/proteinGroups.txt'

#######################################
# Schilling Group data set
#file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2017/11/PXD006914/proteinGroups_normal.txt'
#file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2017/11/PXD006914/proteinGroups_xman.txt'

# nicht so spannend! NRI mit 70.2% only
# file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2017/05/PXD006302/proteinGroups.txt'

# no RI feature
#file = '~/zbsa/2017_RobustQuantilenorm/PublicData/pride/proteinGroups/ftp.pride.ebi.ac.uk/2016/09/PXD003082/txt/proteinGroups.txt'

file = '/Users/schad/zbsa/Schilling/Label-freewithout-EcoliONLY-MBR/txt/proteinGroups.txt'
#######################

str <- unlist(strsplit(file,"/"))
filename <- str[pmatch("PXD",str)]

dat <- read.csv(file, allowEscapes = T, check.names = F,sep = "\t")
expMatrix <- as.matrix(dat[, grepl("^LFQ", names(dat))])
expMatrix[expMatrix == 0] <- NA

#colnames(expMatrix) <- 1:ncol(expMatrix) # A. Schad

allColNA <- as.vector(apply(expMatrix, 1, function(r) {
return(sum(!is.na(r)) == 0)
}))

row.names(expMatrix) <- dat[, "Protein IDs"]
featureAnnotations <- data.frame(proteinName = dat[, "Protein IDs"],
                               proteinDescription = dat[, "Fasta headers"],
                               #idScore = dat[, "Q-value"],
                               isDecoy = dat[, "Reverse"] == "+",
                               nbPeptides = dat[, "Peptides"],
                               isNormAnchor = rep(T, nrow(expMatrix)),
                               isFiltered = dat[, "Reverse"] == "+",
                               row.names = dat[, "Protein IDs"])
##isPotential.contaminant = dat[,"Potential contaminant"]=="+",
featureAnnotations <- featureAnnotations[!allColNA, ]
#expMatrix <- as.matrix(expMatrix[!allColNA, rownames(expDesign)])
expMatrix <- as.matrix(expMatrix[!allColNA, ]) # mod. by A. Schad
expMatrix <- log2(expMatrix)

library("MBQN")
devtools::document()
res <- mbqn.check_saturation(expMatrix, FUN = median, flag_show_fig = TRUE, low_thr = 0.2, filename = filename, feature_index = 124)

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

# use mbqn.boxplot function
# mbqn.boxplot(expMatrix[-452,],nri_max)

# for Schilling data set:
# search for all ECOLI features, show whether their concentration is increasing
#bla <- (expMatrix[grep("COLI",rownames(expMatrix)),])
range(bla,na.rm =T)
#rownames(bla) <- c()
#par(mfrow=c(1,1))
#matplot(t(bla), type = "l")
# now search for Human proteins
# bla <- (expMatrix[grep("HUMAN",rownames(expMatrix)),])
range(bla,na.rm =T)
#rownames(bla) <- c()
#par(mfrow=c(1,1))
#matplot(t(bla), type = "l")
