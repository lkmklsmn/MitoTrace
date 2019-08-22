# testing the github platform

# plot the germline mutation from glioblastoma dataset
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