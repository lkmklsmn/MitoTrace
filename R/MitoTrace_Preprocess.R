# set the working directory
setwd("/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/01_Mito_visualize/02_IGV/")

bams <- list.files("/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/01_Mito_visualize/02_IGV/", full.names = T, pattern = ".bam$")

fasta_loc <- "/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/GitHub/MitoTrace/Data/GRCH38_MT.fa"

# For human
source("/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/GitHub/MitoTrace/R/MitoTrace.R")

SRS1800264 <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")

mae_res <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT")

# test the quality
SRS1800264_q1 <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT", min_base_quality= 1)
SRS1800264_q25 <- MitoTrace(bam_list = bams, fasta = fasta_loc, chr_name = "MT", min_base_quality= 25)

# for plot
mt_ann <- read.csv("/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/18_plot_function/MT_annotation_table.csv", header = TRUE)

MitoPlot(mae=mae_res, species = "human", mt_ann = mt_ann)
