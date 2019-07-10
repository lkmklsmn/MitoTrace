setwd("/Users/mwang14/Documents/01_Zhao/12_Develop_Mitogenotyping_Tools/02_GBM_data")

# Process the GBM scRNA-seq .rds file
raw <- readRDS("processed.MAE_mito.rds")
str(raw)
covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]

af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]

rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], data.frame(rowRanges(allSE))[,c(6)], ">", data.frame(rowRanges(allSE))[,c(7)])

covVARmean <- Matrix::rowMeans(assays(covSE)[["coverage"]][start(rowRanges(allSE)),])
covVARmedian <- matrixStats::rowMedians(data.matrix(assays(covSE)[["coverage"]][start(rowRanges(allSE)),]))
covCell <- Matrix::colMeans(assays(covSE)[["coverage"]])

minCellCoverage <- 30
af2 <- af[covVARmean> 10  & covVARmedian > 10 , covCell > minCellCoverage]
cov2 <- cov[covVARmean> 10 & covVARmedian > 10, covCell > minCellCoverage]
dim(af2)

af3 <- af2[Matrix::rowMeans(af2) > 0.0001 & Matrix::rowMeans(af2) < 0.75, ]
cov3 <- cov2[Matrix::rowMeans(af2) > 0.0001 & Matrix::rowMeans(af2) < 0.75, ]
dim(af3)

# read the SRR patient listed table
ssr_list <- read.table("SRR_Patient_sorted_list",sep = "\t")
head(ssr_list)
dim(ssr_list)
colnames(af3)

length(colnames(af3))
match(colnames(af3), ssr_list$V1)
ssr_list$V3 <- paste(ssr_list$V1, ssr_list$V2, sep = "_")
head(ssr_list)
match(colnames(af3), ssr_list$V1)
ssr_list[(match(colnames(af3), ssr_list$V1)),]
head(ssr_list[(match(colnames(af3), ssr_list$V1)),])
ssr_sorted <- ssr_list[(match(colnames(af3), ssr_list$V1)),]
dim(ssr_sorted)
head(ssr_sorted)
colnames(af3) <-  ssr_sorted$V3
colnames(af3)
dim(af3)
af3[1:5,1:5]

# For the sample name is: BT_S2, BT_S1, BT_S4, BT_S6
#For the first sample
grep("BT_S2", colnames(af3))
x <- NA
for (i in 1:39176){
    x[i] <- t.test(af3[i,1:99], af3[i,100:392])$p.value
}
x <- data.frame(x)
x1 <- cbind(x,rownames(x))
x1_order <- x1[order(x1[,1]),]
attr(x1_order,"names") <- c("pvalue", "row_number")
head(x1_order)
pheatmap::pheatmap(af3[as.numeric(as.character(x1_order$row_number[1:9])),])

pheatmap::pheatmap(af3[as.numeric(as.character(x1_order$row_number[1:9])),], cluster_rows = FALSE, cluster_cols=FALSE)
# 4503  18236 19775 12133 2638  7063  28326 38926 27542

pheatmap::pheatmap(af3[sample(rownames(af3),100),], cluster_rows = FALSE, cluster_cols=FALSE)
sample(rownames(af3),100)

# for the second sample
grep("BT_S1", colnames(af3))
z <- NA
for (i in 1:39176){
  z[i] <- t.test(af3[i,100:199], af3[i,c(1:99,200:392)])$p.value
}
z <- data.frame(z)
z1 <- cbind(z,rownames(z))
z1_order <- z1[order(z1[,1]),]
attr(z1_order,"names") <- c("pvalue", "row_number")
head(z1_order)
pheatmap::pheatmap(af3[as.numeric(as.character(z1_order$row_number[1:10])),])
pheatmap::pheatmap(af3[as.numeric(as.character(z1_order$row_number[1:25])),], cluster_rows = FALSE, cluster_cols=FALSE)
# 25# relative mutation # 15975

# for the third sample
grep("BT_S4", colnames(af3))
k <- NA
for (i in 1:39176){
  k[i] <- t.test(af3[i,200:293], af3[i,c(1:199,294:392)])$p.value
}
k <- data.frame(k)
k1 <- cbind(k,rownames(k))
k1_order <- k1[order(k1[,1]),]
attr(k1_order,"names") <- c("pvalue", "row_number")
head(k1_order)
pheatmap::pheatmap(af3[as.numeric(as.character(k1_order$row_number[1:10])),])

pheatmap::pheatmap(af3[as.numeric(as.character(k1_order$row_number[c(1:3,5:14)])),], cluster_rows = FALSE, cluster_cols=FALSE)
# 30047  2810 23878  15366  9154 14807 26866 31935 34287 21047 8720 36935 8968


# For the last sample
grep("BT_S6", colnames(af3))
y <- NA
for (i in 1:39176){
  y[i] <- t.test(af3[i,294:392], af3[i,1:293])$p.value
}
y <- data.frame(y)
y1 <- cbind(y,rownames(y))
y1_order <- y1[order(y1[,1]),]
head(y1_order)
attr(y1_order,"names") <- c("pvalue", "row_number")
pheatmap::pheatmap(af3[as.numeric(as.character(y1_order$row_number[1:12])),],cluster_rows = FALSE, cluster_cols=FALSE)
# 45  4005  33339 15075 5005  28484 20266 27093

target <- c(4503, 18236, 19775, 12133, 2638, 7063, 28326, 38926, 27542, 15975, 30047, 2810, 23878, 15366, 9154, 14807,26866, 31935, 34287, 21047, 8720, 36935, 8968, 45, 4005, 33339, 15075, 5005, 28484, 20266, 27093)

pheatmap::pheatmap(af3[target,])
pheatmap::pheatmap(af3[target,], cluster_rows = FALSE, cluster_cols=FALSE)

dev.off()
pdf(file = "01_heatmap_BGM.pdf", width = 15, height = 7, paper = "special")
lable <- c(rep("BT_S2",99), rep("BT_S1",100), rep("BT_S4",94), rep("BT_S6",99))
col_fun = colorRamp2(c(0, 1), c("white", "red4"))
ha = HeatmapAnnotation(patient = lable, show_legend = TRUE, col = list(patient = c("BT_S2" = "dodgerblue2", "BT_S1" = "gold1", "BT_S4" = "coral1", "BT_S6" = "darkorchid3")))

Heatmap(target_matrix, 
        col = col_fun, 
        cluster_columns = FALSE, 
        cluster_rows = FALSE, 
        show_column_names = FALSE, 
        border = TRUE, 
        column_title="Mitochondrial germline mutations in human glioblastoma", 
        row_title = "Mitochondiral mutations",
        top_annotation = ha,
        show_heatmap_legend = FALSE
)

lgd = Legend(col_fun = col_fun, title = "Allele\nFrequency", direction = "horizontal", title_position = "topcenter")
draw(lgd, x = unit(2, "cm"), y = unit(0.8, "cm"),just = c("left", "bottom"))
dev.off()


# plot the coverage depth of selected mutations
library("ggplot2")
gbm_depth <- read.table("GBM_coverage_mean_depth",sep = "\t")
site <- c(7775,  827,   2873,  4947,  4924,  11914, 15038, 15535, 13942, 11864, 2092,  5178,  8701,  10873, 15301, 9540,  13105, 4883,  8414,  4790, 14569, 12705, 15043, 709,   6962,  6960,  9950,  8584,  15235, 3537, 13395)
match(site, gbm_depth$V1)
summary(gbm_depth[site,]$V2)
target_depth <- gbm_depth[site,]
theme_set(theme_bw())
ggplot(data=target_depth, aes(x=V1, y=V2)) + 
  geom_point(col="blue") + 
  geom_line(col="blue") + 
  ggtitle("Selected MT germline mutation mean depth")+
  theme(plot.title = element_text(hjust = 0.5))+
  ylab("Mean_coverage_depth")+
  xlab("MT genome position")+
  ylim(0,20)

      
# plot the MT germline mutation frequency
target_mutation <- c("G>A","A>G","A>G","T>C","G>A","G>A","A>G","C>T","A>G","T>C","C>T","C>A","A>G","T>C","G>A","T>C","A>G","C>T","C>T","A>G","G>A","C>T","G>A","G>A","G>A","C>T","T>C","G>A","A>G","A>G","A>G")
table(target_mutation)
barplot(table(target_mutation))
barplot(table(target_mutation), anlge = c(45), density = 20, col = "blue", ylim=c(0,12), main = "MT germline mutation frequency")