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

The package provides the basic function: `mbqn()` which does quantile normalization or mean balance quantile normalization of a matrix object. The matrix may contain NAs. The argument `FUN` is used to select between classical quantile normalization, and mean or median balanced quantile normalization.

## Examples
1. Generate simple matrix, apply median-balanced quantile normalization, generate boxplot of normalized features and check for rank invariant (RI) or nearly rank invariant (NRI) features:

`mtx <-  matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3)`<br/>
`mtx.norm <- mbqn(x = mtx, FUN = median)`<br/>
`mbqn.boxplot(mtx.norm)`<br/>
`mbqn.check_saturation(mtx,FUN = median,low_thr = 0.1, filename = "simple_mtx",feature_index = 1)`

2. This example will download data from the PRIDE repository, normalize the data, identifies RI/NRI features, and give graphical output. The example is found in folder /installationpath/mbqn/examples/.

`example1(3)`

To run example1.R for data from the PRIDE archive, the respective proteinGroups.txt file must be downloaded from the PRIDE webpage to to the folder /installationpath/mbqn/examples/PXDxxxx/ or directly by the function in exmple1.R which uses the R package rpc. 

## Figures
Figures created by MBQN are saved under /installationpath/mbqn/.
