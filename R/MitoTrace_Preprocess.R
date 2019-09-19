# set the working directory
setwd("/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/03_develop_script/01_pileup_Rsamtools/")

bams <- list.files("/Users/mwang14/Documents/01_Zhao/01_raw_bam/", full.names = T, pattern = ".bam$")
fasta_loc <- "/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/03_develop_script/01_pileup_Rsamtools/GRCH38_MT.fa"
