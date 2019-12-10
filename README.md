# MitoTrace
This repository contains the R code for MitoTrace, a computational framework to study mitochondrial variation from single-cell RNA sequencing data.


## Prerequisites
MitoTrace runs on 32-bit or 64-bit GNU/Linux R environment and requires the following dependencies: 
* R (>= 3.6.1)
* seqinr (>= 3.4-5)
* Matrix (>= 1.2-17)
* Rsamtools (>= 2.0.0)

Please  make sure these packages (and correct versions) are installed or install yourself in the following way:
```
install.packages("seqinr")
install.packages("Matrix")
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsamtools")
```


## Install
If you want to use the lastest version of MitoTrace, please install it directly from GitHub. First, you need to install and load the devtools package. Then use install_github("lkmklsmn/MitoTrace") as follow:
```
install.packages("devtools")
library("devtools")
install_github("lkmklsmn/MitoTrace")
```

## Usage

Plot the read coverage across the mitochondria
```
MitoDepth(bam_list = bams, species = "human", mt_ann = mt_ann)
```

Extract read coverage and alternative allele counts across all positions
```
mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")
```

Calculate allele frequency
```
af <- calc_allele_frequency(mae_res)
```
## Examples
Please check the *examples* folder for Markdown files, or use the following link to open.

1. [Reproduce the Leif et al analysis[1]](http://htmlpreview.github.com/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Reproduce_Cell_Leif_et_al.html)

2. [Identify personal variants in SMART-seq2 data[2]](http://htmlpreview.github.com/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Single-Cell-SMART-SEQ2-data.html)

3. [Distinguish different cell lines using 10X Genomics data](http://htmlpreview.github.com/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Single-Cell-10X-Genomics-data.html)

## Citation

## License
`MitoTrace` uses GNU General Public License GPL-3.

## Reference
1.	Ludwig, L.S., et al., Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics. Cell, 2019. 176(6): p. 1325-1339 e22.
2.	Darmanis, S., et al., Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. Cell Rep, 2017. 21(5): p. 1399-1410.

### Note
Our opimized, efficient and user-friendly R package `MitoTrace` is currently under development.
