## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ----gh-installation, eval = FALSE---------------------------------------
#  # install.packages("devtools")
#  devtools::install_github("arianeschad/mbqn")

## ----installation, eval = FALSE------------------------------------------
#  # install.packages("githubinstall")
#  githubinstall::githubinstall("mbqn")

## ----bc-installation, eval = FALSE---------------------------------------
#  # if (!requireNamespace("BiocManager", quietly = TRUE))
#  #    install.packages("BiocManager")
#  BiocManager::install("mbqn")

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install("limma")

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install("preprocessCore")

## ---- eval = FALSE-------------------------------------------------------
#  BiocManager::install("rpx")

## ---- eval = FALSE-------------------------------------------------------
#  install.packages(pkgs = c("filesstrings"), dependencies = TRUE)

## ----example1, eval = TRUE-----------------------------------------------
## basic example
library("MBQN")
# mtx <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3)
mtx <- mbqnSimuData("omics.dep")
mtx <- mbqnSimuDistortion(mtx)$x.mod

## ----figure1, fig.height = 5, fig.width = 6, fig.align = "center", fig.cap = "Fig. 1 Boxplot of the unnormalized, distorted intensity data matrix. The first feature is an RI feature (red line). It has maximum intensity for each sample!"----
plot.new()
mbqnBoxplot(mtx, irow = 1, main = "unnormalized")

## ---- eval = TRUE--------------------------------------------------------
ri.obj <- mbqnGetNRIfeatures(mtx, low_thr = 0.5)

## ----figure2, fig.height = 3, fig.width = 6, fig.align = "center", fig.cap = "Fig. 2 Correctly detected and identified RI feature with a data coverage of 100% across samples."----
plot.new()
mbqnPlotRI(ri.obj)

## ----figure3, fig.height = 5, fig.width = 6, fig.align = 'center', fig.cap = "Fig. 3 Quantile normalized intensities with balanced and unbalanced normalized RI feature. Classical quantile normalization suppresses any intensity variation of the RI feature, while the MBQN preserves its variation while reducing systematic batch effects!"----
plot.new()
mtx.norm <- mbqnNRI(x = mtx, FUN = median, verbose = FALSE) # MBQN
mtx.qn <- mbqnNRI(x = mtx, FUN = NULL, verbose = FALSE) # QN
mbqnBoxplot(mtx.norm, irow = ri.obj$ip, vals = data.frame(QN = mtx.qn[ri.obj$ip,]), main = "normalized")

## ----example2, eval = FALSE----------------------------------------------
#  ## Normalize LFQ intensity data from the PRIDE repository
#  mbqnExample(which.example = 1)

## ---- eval = FALSE-------------------------------------------------------
#  ??MBQN

