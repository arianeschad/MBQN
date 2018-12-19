
<!-- README.md is generated from README.Rmd. Please edit that file -->
MBQN Package
============

Mean/Median-balanced quantile normalization for processing omics data

Description
-----------

This package supplies a modified quantile normalization for preprocessing omics or other matrix-like organized data with intensity values biased by global shifts of mean and scatter of columns. The modification balances the mean intensity of features which are rank invariant (RI) or nearly rank invariant (NRI), e.g. features that have always largest intensity across samples/columns, before quantile normalization in order to reduce systematics in downstream analysis.

Installation
------------

To install this package from Github, you need R version &gt;= 3.3.3.

For installation from source run in R:

``` r
# install.packages("devtools")
devtools::install_github("arianeschad/mbqn")
```

or

``` r
# install.packages("githubinstall")
githubinstall::githubinstall("mbqn")
```

Installation from source via CRAN: <br/>

``` r
install.packages("mbqn")
```

Dependencies:
-------------

The package uses `normalize.quantiles()` from the package preprocessCore by B. Bolstad (2016), preprocessCore: A collection of pre-processing functions. R package version 1.36.0. available from <https://github.com/bmbolstad/preprocessCore>: <br/> `install.packages("preprocessCore")`

Optionally, the limma package can be used for computation of the quantile normalization: <br/>

``` r
install.packages("limma")
```

Collecting data from PRIDE experiments in `example1()` requires the R package rpx by L. Gatto (2017), version 1.10.2. available from <https://github.com/lgatto/rpx>: <br/>

``` r
install.packages("rpx")
```

In R run <br/>

``` r
install.packages(pkgs = c("preprocessCore","limma","rpx"), dependencies = TRUE)
```

Additional packages needed to run MBQN examples: <br/>

``` r
install.packages(pkgs = c("filesstrings"), dependencies = TRUE)
```

After installation, check the `?mbqn.demo` and the `?mbqn.example1` help for full working examples with further documentation.

Basic Usage
-----------

The package provides two basic functions: `mbqn()` applies quantile normalization or mean-balanced quantile normalization to a matrix. `mbqn.nri()` applies quantile normalization and mean balanced quantile normalization only to selected nearly rank invariant and rank invariant features, specified by a threshold or manually. The input matrix may contain NAs. To run one of these functions you will need to provide an input matrix similar to the data matrix in `mbqn.demo` or `mbqn.example1`. The argument `FUN` is used to select between classical quantile normalization (default), and mean or median balanced quantile normalization. The function `mbqn.check_saturation()` can be used to check a data matrix for rank or nearly rank invariant features. It provides a list of potential RI/NRI features, a rank invariance frequency, and a graphical output.

Examples
--------

Example 1: Generate a simple matrix, apply median-balanced quantile normalization, generate a boxplot of normalized features and check for rank invariant (RI) or nearly rank invariant (NRI) features:

``` r
## basic example
mtx <- matrix(c(5,2,3,NA,4,1,4,2,3,4,6,NA),ncol=3)
mtx.norm <- MBQN::mbqn(x = mtx, FUN = median)
MBQN::mbqn.boxplot(mtx.norm, irow = 1,)
MBQN::mbqn.check_saturation(mtx,FUN = median,low_thr = 0.1, filename = "simple_mtx",feature_index = 1)
```

Example 2: This example will download data with LFQ intensities from the PRIDE repository, normalize the data, identifies RI/NRI features, and give graphical output. One can choose between four data sets. By default data files are stored in the current working directory currentdir/examples/PXDxxx.

``` r
## Normalize LFQ intensity data from the PRIDE repository
mbqn.example1(which.example = 3)
```

In order to run `mbqn.example1()`, the respective proteinGroups.txt file must be first downloaded from the PRIDE repository to the folder currentdir/examples/PXDxxxx/ in advance or you can directly download the data by running the code included in `mbqn.example1()`.

Figures
-------

Figures created with MBQN are saved in the current working directory.

References
----------

A. Schad and C. Kreutz, MBQN: R package for mean balanced quantile normalization. In prep. for Bioinf. Appl. Note, 2018
