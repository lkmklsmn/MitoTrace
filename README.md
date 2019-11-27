# MitoTrace
This repository contains the R code for MitoTrace. A R package to infer mitochondrial heterplasmies from single-cell RNA sequencing data. 


## Dependencies for MitoTrace
* R (>= 3.6.0)
* seqinr (>= 3.4-5)
* Matrix (>= 1.2-17)
* Rsamtools (>= 2.0.0)


## Install

`MitoTrace` is available on Bioconductor, you could install it by:

```
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("MitoTrace")
```
If you want to use the lastest version, please install it directly from GitHub:
```
library("devtools")
install_github("lkmklsmn/MitoTrace")
```


## Usage


## Example

##### `MitoTrace::DepthPlot` draw the coverage depth of imported single cell data
![GitHub Logo](https://github.com/lkmklsmn/MitoTrace/blob/master/example/gene_bar_cov.png)

##### `MitoTrace` enables identification of germline variants from single cell SMART-SEQ2 data
![GitHub Logo](https://github.com/lkmklsmn/MitoTrace/blob/master/example/R1.png)

##### `MitoTrace` enables identification of germline variants from single cell 10X genomics data
![GitHub Logo](https://github.com/lkmklsmn/MitoTrace/blob/master/example/10x_genomics.png)

```
# read the file
demuxlet <- fread("jurkat_293t_demuxlet.best")
bams <- list.files("/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/21_read_10xGenomics/jurkat", full.names = T, pattern = ".bam$")
fasta_loc <- "/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/GitHub/MitoTrace/Data/GRCH38_MT.fa"

mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")
  
# Run MitoTrace
mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT", min_umi = 100)

# Calculate allele frequencies
af <- calc_allele_frequency(mae_res)

# Generate plots
ok <- intersect(colnames(af), demuxlet$BARCODE)
af <- af[,ok]
cell_line <- sapply(demuxlet$BEST,function(x){strsplit(x,"-")[[1]][[2]]})
cell_line <- cell_line[match(ok, demuxlet$BARCODE)]

good_variants <- names(which(rowMeans(af) > 0.3))
good_variants <- names(tail(sort(apply(af, 1, var)), 20))
af_hv <- data.matrix(af[good_variants,])
ydata <- Rtsne::Rtsne(t(af_hv))

aframe <- data.frame(cell_line, ydata$Y)
ggplot(aes(X1, X2, color = cell_line), data = aframe) + geom_point()

anno <- data.frame(cell_line)
rownames(anno) <- colnames(af_hv)
pheatmap::pheatmap(af_hv, show_colnames = F, annotation_col = anno)
```


##### `MitoTrace` detects increased mutational burden
![GitHub Logo](https://github.com/lkmklsmn/MitoTrace/blob/master/example/barplot.png)

## Citation


## License
`MitoTrace` uses GNU General Public License GPL-3.

## Reference
1.	Ludwig, L.S., et al., Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics. Cell, 2019. 176(6): p. 1325-1339 e22.
2.	Darmanis, S., et al., Single-Cell RNA-Seq Analysis of Infiltrating Neoplastic Cells at the Migrating Front of Human Glioblastoma. Cell Rep, 2017. 21(5): p. 1399-1410.
3.	Angelidis, I., et al., An atlas of the aging lung mapped by single cell transcriptomics and deep tissue proteomics. Nat Commun, 2019. 10(1): p. 963.

### Note
Our opimized, efficient and user-friendly R package `MitoTrace` is currently under development.