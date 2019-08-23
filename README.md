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

##### `MitoTrace` enables identification of germline variants

##### `MitoTrace` enables identification of somatic variants

##### `MitoTrace` detects increased mutational burden


## Citation


## License
`MitoTrace` uses GNU General Public License GPL-3.

## Reference
1.	Ludwig, L.S., et al., Lineage Tracing in Humans Enabled by Mitochondrial Mutations and Single-Cell Genomics. Cell, 2019. 176(6): p. 1325-1339 e22.

### Note
Our opimized, efficient and user-friendly R package `MitoTrace` is currently under development.