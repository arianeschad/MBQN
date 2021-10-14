
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MBQN Package

Mean/Median-balanced quantile normalization for preprocessing omics data

## Description

This package contains a modified quantile normalization (QN) for
preprocessing and analysis of omics or other matrix-like organized data
with intensity values biased by global, columnwise distortions of
intensity mean and scale. The modification balances the mean intensity
of features (rows) which are rank invariant (RI) or nearly rank
invariant (NRI) across samples (columns) before quantile normalization
\[1\]. This helps to prevent an over-correction of the intensity
profiles of RI and NRI features by classical QN and therefore supports
the reduction of systematics in downstream analyses. Additional package
functions help to detect, identify, and visualize potential RI or NRI
features in the data and demonstrate the use of the modification.

## Installation

To install this package, you need R version \>= 3.6.

For installation from Bioconductor run in R:

``` r
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
BiocManager::install("MBQN")
```

For installation from Github run in R:

``` r
install.packages("devtools")
devtools::install_github("arianeschad/MBQN")
```

or

``` r
install.packages("githubinstall")
githubinstall::githubinstall("MBQN")
```

## Dependencies

The core of the MBQN package uses `normalizeQuantiles()` from the
package `limma`\[2\], available at
<https://bioconductor.org/packages/release/bioc/html/limma.html>, for
computation of the quantile normalization. Optionally,
`normalize.quantiles()` from the package `preprocessCore`\[3\],
available at
<https://bioconductor.org/packages/release/bioc/html/preprocessCore.html>,
can be used. <br/>

The function `getPXDfile()` in MBQN uses data from the PRIDE repository.
To run this function one needs the R package `rpx` \[4\] to download the
data. <br/>

To install these packages in R run: <br/>

``` r
# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# BiocManager::install(pkgs = c("preprocessCore","limma","rpx","SummarizedExperiment"))
```

## Usage

Further information about the package is provided at the wiki <br/>
<https://github.com/arianeschad/MBQN/wiki>

## References

\[1\] Brombacher, E., Schad, A., Kreutz, C. (2020). Tail-Robust Quantile
Normalization. Proteomics. <br/> \[2\] Ritchie, M.E., Phipson, B., Wu, D.,
Hu, Y., Law, C.W., Shi, W., and Smyth, G.K. (2015). limma powers
differential expression analyses for RNA-sequencing and microarray
studies. Nucleic Acids Research 43(7), e47. <br/> \[3\] Ben Bolstad
(2018). preprocessCore: A collection of pre-processing functions. R
package version 1.44.0. <https://github.com/bmbolstad/preprocessCore>.
<br/> \[4\] Laurent Gatto (2019). rpx: R Interface to the
ProteomeXchange Repository. R package version 1.18.1.
<https://github.com/lgatto/rpx>.
