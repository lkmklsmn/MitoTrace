# MitoTrace
Mitochondrial variation is linked to important biological functions and various human diseases. Recent progress in single-cell genomics has established single-cell RNA sequencing (scRNAseq) as a popular and powerful technique to profile transcriptomics at the cellular level. Here, we developed MitoTrace, an R package for the analysis of mitochondrial variation in scRNAseq data. This repository contains the R code for MitoTrace.

## Prerequisites

MitoTrace runs on 32-bit or 64-bit GNU/Linux R environment and requires the following dependencies: \* R (\>= 3.6.1) \* seqinr (\>= 3.4-5) \* Matrix (\>= 1.2-17) \* Rsamtools (\>= 2.0.0)

Please make sure these packages (and correct versions) are installed or install yourself in the following way:

    install.packages("seqinr")
    install.packages("Matrix")
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("Rsamtools")

## Install
MitoTrace is designed for the R programming language  and statistical computing environment. If you want to use the lastest version of MitoTrace, please install it directly from GitHub. First, you need to install and load the devtools package. Then use install_github("lkmklsmn/MitoTrace") as follow line in your R console:

    install.packages("devtools")
    library("devtools")
    install_github("lkmklsmn/MitoTrace")

## Usage

MitoTrace plots the read coverage across the mitochondria human or mouse genome

    MitoDepth(bam_list = bams, species = "human", mt_ann = mt_ann)

MitoTrace extracts read coverage and alternative allele counts across all mitochondrial genome positions

    mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")

MitoTrace calculates allele frequency for each mitochondrial position

    af <- calc_allele_frequency(mae_res)

## Examples

Please check the *examples* folder for Markdown files, or use the following link to open.

1.  [MitoTrace reproduces the heteroplamsy reported in Leif's lineage tracing study analysis[1]](https://htmlpreview.github.io/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Reproduce%20results%20in%20figures%205H%20%26%205I%20from%20Leif%20et%20al.html)

2.  [MitoTrace identifies personal variants in SMART-seq2 data data[2]](https://htmlpreview.github.io/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Single-Cell-SMART-SEQ2-data.html)

3.  [MitoTrace distinguishes different cell lines in 10X Genomics data data[3]](https://htmlpreview.github.io/?https://github.com/lkmklsmn/MitoTrace/blob/master/examples/Single-Cell-10X-Genomics-data.html)

## License

`MitoTrace` uses GNU General Public License GPL-3.

## References

1.  Ludwig, L.S., et al., Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics. Cell, 2019.
2.  Darmanis, S., et al., Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. Cell Rep, 2017. 
3.  McGinnis, G., et al., DoubletFinder: Doublet Detection in Single-Cell RNA Sequencing Data Using Artificial Nearest Neighbors. Cell Systems, 2019.

### Note

Our optimized, efficient and user-friendly R package `MitoTrace` is currently under development.
