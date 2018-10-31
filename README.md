---
title: "MBQN Package"
output: github_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```
# mbqn
Mean/Median-balanced quantile normalization for processing omics data

## GitHub Documents

This packages provides a modified quantile normalization function for omics data or other matrix-like data. The normalization removes systematic batch effects between measurements. The modification is consists of a mean balancing which reduces systematics in downstream analysis for features that are always or mostly of largest intensity accross all samples. This function uses normalize.quantiles() from the package preprocessCore that can be installed from http://bioconductor.org/biocLite.R
When you click the **Knit** button all R code chunks are run and a markdown file (.md) suitable for publishing to GitHub is generated.

### Installing The Package

To install this package from Github, you will need to Hadley Wickham's devtools package installed.

install.packages("devtools")
library("devtools")

Now we can install from Github using the following line:

devtools::install_github("matthewjdenny/ContentStructure")

I have  had success installing with R 3.2.0+ installed but if you do not have the latest version of R installed, or run into some install errors (please email if you do), it should work as long as you install the dependencies first with the following block of code:

install.packages( pkgs = c("BH","coda","RcppArmadillo","gridBase",
"gplots","slam","vegan"), dependencies = TRUE)

If all went well, check out the `?ContentStructure` help file to see a full working example with info on how the data should look.

## Basic Useage

The package provides the basic function: `mbqn()` which does the normalization. To run this function, you will need to provide input similar to the example `example1` included as examples with the package.

## Example

This example code will download data from the PRIDE repository, normalize the data, identifies RI/NRI features and give graphical output. 

# Load in necessary data
example1(3)

Note that the `echo = FALSE` parameter was added to the code toxxx .
