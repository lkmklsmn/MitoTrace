# set directory
setwd("C:/Users/mwang14/Google Drive/00_Texas_Posdoc_Career/01_Zhao/12_Develop_Mitogenotyping_Tools/18_plot_function/01_ori_coverage_depth/")
human_srr <- read.table("SRR_human_mean_depth", header = FALSE)


# prepare the color
library("RColorBrewer")
display.brewer.all()
display.brewer.pal(n = 12, name = 'Set3')
col= c(brewer.pal(n = 12, name = "Set3"), brewer.pal(n = 11, name = "Spectral"), brewer.pal(n = 9, name = "Set1"), brewer.pal(n = 8, name = "Accent"))
col

pdf("04_gene_bar_cov_plot.pdf", width=25, height = 12.5, paper = "special")
plot(human_srr$V1, human_srr$V2, type = "l", col="brown1", xlab="MT genome postion", ylab="log2(Coverage depth)",main="scRNA-seq coverage depth cross MT genome", ylim=c(-5,21))

# set the gene label location
location <- c(288.55, 577, 1124.55, 1602, 2450.05, 3230, 3784.55, 
              3900, 4300, 4700, 
              4990.55, 
              5150, 5550, 5950, 6350, 6750, 7150.55, 
              7300, 7718, 
              7927.55, 8295, 8469.05, 8867.05, 9598.55,
              9991, 10231.55, 10405, 10618.05, 11448.55, 
              11750, 12200, 12650, 
              13242.55, 14411.05, 14674, 15317.05, 
              15700, 16200, 16584.55)

# using rect function instead of segments
for (i in 1:39) {
  if(grepl("TR", mt_ann$gene[i])){
    rect(mt_ann[i,2], -2.5, mt_ann[i,3], -1, col=col[i], border = "white")
    text(location[i], -0.5, mt_ann$gene[i], cex = 0.9)
  }
  else if(grepl("RN", mt_ann$gene[i])){
    rect(mt_ann[i,2], -2.5, mt_ann[i,3], -1, col=col[i], border = "white")
    text(location[i], -0.5, mt_ann$gene[i], cex = 0.9)
  }
  else{
    rect(mt_ann[i,2], -4, mt_ann[i,3], -2.5, col=col[i], border = "white")
    text(location[i], -4.5, mt_ann$gene[i], cex = 0.9)
  }
}

dev.off()
