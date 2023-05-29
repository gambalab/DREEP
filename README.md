# Drug Response Estimation from single-cell Expression Profiles (DREEP)

## Description

DREEP (Drug Response Estimation from single-cell Expression Profiles) is a bioinformatics tool that leverages results from large-scale cell-line viability screens and enrichment analysis to predict drug vulnerability from the transcriptional state of a cell. It only requires a pre-defined collection of Genomic Profiles of Drug Sensitivity (GPDS) signatures that are ranked lists of genes reflecting their importance in predicting the effect of a small molecule.

# Installation

DREEP require the installation of [gficf v2](https://github.com/gambalab/gficf) package first. Follow the instruction below to install it on Ubuntu. For other OS althoug not officially supported you can give a look [here](https://htmlpreview.github.io/?https://github.com/gambalab/gficf/blob/master/inst/doc/installation.html)

### 1. Installation of GFICF on Ubuntu/Debian

You need gsl dev library to successfully install RcppGSL library. On Ubuntu/Debian systems this can be accomplished by running the following command `sudo apt-get install libgsl-dev libcurl4-openssl-dev libssl-dev libxml2-dev` from the terminal.

Then exec in R terminal the following commands to install gficf

``` r
# Install required bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {install.packages("BiocManager")}

BiocManager::install(setdiff(c("sva","edgeR", "fgsea"),rownames(installed.packages())),update = F)

# We rquire RcppML package from github (not the cran version)
if("RcppML" %in% rownames(installed.packages())) {remove.packages("RcppML")}
devtools::install_github("zdebruine/RcppML")

if(!require(devtools)){ install.packages("devtools")}
devtools::install_github("gambalab/gficf")
```

### 2. Installation of DREEP

From the R terminal the following commands

``` r
if(!require(devtools)){ install.packages("devtools")}
devtools::install_github("gambalab/DREEP")
```

# Example of use



