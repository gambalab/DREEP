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

From the R terminal exec the following commands

``` r
if(!require(devtools)){ install.packages("devtools")}
devtools::install_github("gambalab/DREEP")
```

### 3. Example of use
``` r
library(gficf)
library(DREEP)
library(ggplot2)

data(small_BC_atlas)
data <- gficf(M=small_BC_atlas,verbose = T)

# Run DREEP on all the cell of the atlas using only CTRP2 and GDSC drug datasets
dereep.data <- DREEP::runDREEP(M = data$gficf,n.markers = 250,gsea = "multilevel",gpds.signatures = c("CTRP2","GDSC"))

# DREEP predictions are into dereep.data$df data frame
# Each row is a drug and the coloumn sens contains
# the percentage of cells predicted sensitive to the drug
head(dereep.data$df)


# Run cell reduction using drug response profile estimated from DREEP
dereep.data <- DREEP::runDrugReduction(dereep.data,verbose = F,cellDistAbsolute = T,reduction = "umap")

# UMAP coordinate are stored into dereep.data$embedding data frame
head(dereep.data$embedding)

# Let's add cell line names and plot the recontructed cell embedding space
dereep.data$embedding$ccl <- sapply(strsplit(x = rownames(dereep.data$es.mtx),split = "_",fixed = T),function(x) x[1])
ggplot(data = dereep.data$embedding,aes(x=X,y=Y,color=ccl)) + geom_point(size=.5) + theme_bw() + xlab("UMAP 1") + ylab("UMAP 2")
```



