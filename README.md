# MitoTrace
This repository contains the R code for MitoTrace. A R package to infer mitochondrial heterplasmies from single-cell RNA sequencing data. 


Install 

MitoTrace is available on Bioconductor, you could install it by:

if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("MitoTrace")

If you want to the lastest version, please install it directly from GitHub:
library("devtools")
install_github("lkmklsmn/MitoTrace")