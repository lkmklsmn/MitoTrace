# MitoTrace
This repository contains the R code for MitoTrace, a computational framework to infer mitochondrial heterplasmies from single-cell RNA sequencing data.


## Prerequisites
MitoTrace runs on 32-bit or 64-bit GNU/Linux R environment and requests dependent R packages: seqinr, Matrix, Rsamtools. Here is are the concrete version of the each dependencies.
* R (>= 3.6.1)
* seqinr (>= 3.4-5)
* Matrix (>= 1.2-17)
* Rsamtools (>= 2.0.0)

Here are the command line used to install the packages from MitoTrace.

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
MitoDepth(bam_list = bams, species = "human", mt_ann = mt_ann)

mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")

af <- calc_allele_frequency(mae_res)

## Examples
Please check the *examples* folder for Markdown files, or use the following link to open.

Reproduce the Leif et al cell result: <http://htmlpreview.github.com/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Reproduce_Cell_Leif_et_al.html>

Single-cell SMART-SEQ2: 
<http://htmlpreview.github.com/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Single-Cell-SMART-SEQ2-data.html> 

Single-cell 10X Genomics: 
<http://htmlpreview.github.com/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Single-Cell-10X-Genomics-data.html>

## Citation

## License
`MitoTrace` uses GNU General Public License GPL-3.

## Reference
1.	Ludwig, L.S., et al., Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics. Cell, 2019. 176(6): p. 1325-1339 e22.
2.	Darmanis, S., et al., Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. Cell Rep, 2017. 21(5): p. 1399-1410.
3.	Angelidis, I., et al., An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics. Nat Commun, 2019. 10(1): p. 963.

### Note
Our opimized, efficient and user-friendly R package `MitoTrace` is currently under development.
