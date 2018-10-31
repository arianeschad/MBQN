# MBQN Package
Mean/Median-balanced quantile normalization for processing omics data

## Description
This package provides a modified quantile normalization function for omics data or other matrix-like data. The modification consists of a mean balancing which reduces systematics in downstream analysis for features that are always or mostly of largest intensity accross all samples. 

## Installing the Package

To install this package from Github, you need R version >= 3.3.3 and the package devtools or githubinstall.

In R use the following line:<br/>
`install.packages("devtools")`<br/>
`devtools::install_github("arianeschad/mbqn")`

or:

`install.packages("githubinstall")`<br/>
`githubinstall::githubinstall("mbqn")`

## Additional dependencies: 
This function uses normalize.quantiles() from the package preprocessCore that can be installed from http://bioconductor.org/biocLite.R It is installed by<br/>
`install.packages("preprocessCore")`

Optionally, the limma package can be used for computation of the quantile normalization. It is installed by <br/>
`install.packages("limma")`

Collecting data from PRIDE experiments in `example1()` requires the rpx package
`install.packages("rpx")`
by Laurent Gatto (2017). rpx: R Interface to the ProteomeXchange Repository. R package version 1.10.2. https://github.com/lgatto/rpx


# Basic Useage

The package provides two basic functions: `mbqn()` applies quantile normalization or mean-balanced quantile normalization to a matrix. The matrix may contain NAs. The argument `FUN` is used to select between classical quantile normalization (default), and mean or median balanced quantile normalization. The function `mbqn.check_saturation()` can be used to check a data matrix for rank or nearly rank invariant features. It provides a list of potential RI/NRI features, a rank invariance frequency, and a graphical output. 

## Examples
1. Generate a simple matrix, apply median-balanced quantile normalization, generate a boxplot of normalized features and check for rank invariant (RI) or nearly rank invariant (NRI) features:

`mtx <-  matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3)`<br/>
`mtx.norm <- mbqn(x = mtx, FUN = median)`<br/>
`mbqn.boxplot(mtx.norm, irow = 1)`<br/>
`mbqn.check_saturation(mtx,FUN = median,low_thr = 0.1, filename = "simple_mtx",feature_index = 1)`

2. This example will download data with LFQ intensities from the PRIDE repository, normalize the data, identifies RI/NRI features, and give graphical output. The example is found in the folder /installationpath/mbqn/examples/.

`example1(which.example = 3)`

To run `example1()` for data from the PRIDE archive, the respective proteinGroups.txt file must be downloaded from the PRIDE webpage to the folder /installationpath/mbqn/examples/PXDxxxx/ in advance or directly by running the code included in `exmple1()` which uses the R package rpc. One can choose between four selected data sets. 

## Figures
Figures created by MBQN are saved under /installationpath/mbqn/.

## References
Please cite: A. Schad and C. Kreutz, MBQN: R package for mean balanced quantile normalization. In prep. for Bioinf. Appl. Note, 2018

