
# Read in MitoTrace output ####
raw <- readRDS('/Users/lukas.simon/OneDrive/Miko/UTHealth/MitoTrace/FromTim/GBM_processed.MAE_mito.rds')
covSE <- raw[["coverage"]]
allSE <- raw[["alleles"]]
af <- assays(allSE)[["counts"]]/(assays(covSE)[["coverage"]][start(rowRanges(allSE)),] + 0.001)
cov <- assays(covSE)[["coverage"]][start(rowRanges(allSE)),]
rownames(af) <- paste0(data.frame(rowRanges(allSE))[,c(2)], data.frame(rowRanges(allSE))[,c(6)], ">", data.frame(rowRanges(allSE))[,c(7)])
info <- read.delim("/Users/lukas.simon/OneDrive/Miko/UTHealth/MitoTrace/DarmanisEtAl/GBM_SRA_plate_relation")
colnames(cov) <- colnames(af) <- info$plate.well[match(colnames(af), info$sra)]

# Read in patient labeling ####
ssr_list <- read.table("/Users/lukas.simon/OneDrive/Miko/UTHealth/MitoTrace/FromTim/SRR_Patient_all_list_sorted",sep = "\t")
ssr_list <- ssr_list[match(colnames(af), ssr_list$V1),]
geo <- read.delim("/Users/lukas.simon/OneDrive/Miko/UTHealth/MitoTrace/DarmanisEtAl/FromGeoSpreadsheet.tsv", row.names = 1)
geo <- geo[colnames(af),]

# Generate tSNE of highly variable SNPs ####
vars <- apply(af, 1, var)
ok <- names(tail(sort(vars), 500))
tmp <- data.matrix(af[ok,])
tmp <- tmp[,which(apply(tmp, 2, var) > 0)]
correl <- cor(tmp)
tData <- Rtsne::Rtsne(1 - abs(correl), perplexity = 20)
plot(tData$Y, col = as.character(meta$Sample.name.color), main = "Patient", pch =20)

# Split by allele frequency ####
asplit <- split(ssr_list$V1, ssr_list$V2)
sample_af_means <- do.call(cbind, lapply(asplit, function(x){
  subm_af <- af[,x]
  #subm_cov <- cov[,x]
  Matrix::rowMeans(subm_af)
}))
cutoff <- 0.75
ok <- which(apply(sample_af_means, 1, function(x) sum(x > cutoff)) > 0)

# Identify informative germline mutations using AOV ####
ok <- which(rowMeans(sample_af_means) > 0.1)
pvals <- apply(af[ok,], 1, function(x) summary(aov(x ~ geo$V8))[[1]][1, 5])
ok <- names(which(pvals < 1e-100))

# Generate heatmap ####
asplit <- split(colnames(af), geo$V8)
samples_ok <- unlist(lapply(asplit, function(x) sample(x, 200)))
tmp <- af[ok, samples_ok]
anno_col <- data.frame(patient = geo$V8[match(samples_ok, rownames(geo))])
rownames(anno_col) <- colnames(tmp)

pheatmap::pheatmap(tmp, show_colnames = F, annotation_col = anno_col, cluster_cols = T)


##########################
ok <- which(geo$V8 == "BT_S4")
geo_ok <- geo[ok,]
asplit <- split(rownames(geo_ok), geo_ok$V10 == "Neoplastic")
sample_af_means <- do.call(cbind, lapply(asplit, function(x){
  subm_af <- data.matrix(af[,x])
  subm_cov <- data.matrix(cov[,x])
  subm_af[which(subm_cov < 20)] <- NA
  rowMeans(subm_af[, sample(colnames(subm_af), 100)], na.rm = T)
}))
tmp <- af[,unlist(asplit)]
ok <- which(apply(tmp, 1, function(x) sum(x > 0.1)) > 1)
treat <- geo[colnames(tmp), "V10"] == "Neoplastic"
pvals <- apply(tmp[ok,], 1, function(x) summary(aov(x ~ treat))[[1]][1,5])

boxplot(split(tmp["1412G>A",], treat))
ok <- which(abs(delta) > 0.25)
#d <- dist(t(cnv[ok, colnames(tmp)]))
#h <- hclust(dist(d))

meta <- meta[colnames(tmp),]
pca <- svd(cnv[sample(rownames(cnv), 5000), colnames(tmp)])
library(mclust)
plot(pca$v, col = as.character(meta$Cluster_2d_color), xlab = "PC1", ylab = "PC2")
legend("bottomright", c("Neoplastic", "Non-neoplastic"), col = c("green", "red"))

plot(pca$v, col = (tmp["1412G>A",] > 0) + 1, xlab = "PC1", ylab = "PC2")

clustering <- Mclust(pca$v[,1:2], G = 2)
plot(pca$v, col = (clustering$classification) + 1, xlab = "PC1", ylab = "PC2")
ok <- which(apply(tmp, 1, function(x) sum(x > 0.1)) > 1)
treat <- clustering$classification
pvals <- apply(tmp[ok,], 1, function(x) summary(aov(x ~ treat))[[1]][1,5])

plot(pca$v, col = (tmp["1412G>A",] > 0) + 1, xlab = "PC1", ylab = "PC2")

vars <- apply(tmp, 1, var)
r <- names(tail(sort(vars), 200))
pca2 <- svd(tmp[r,])
plot(pca2$v, col = (clustering$classification) + 1, xlab = "PC1", ylab = "PC2")
clustering2 <- kmeans(pca2$v[,1:2], centers = 2)

pheatmap::pheatmap(tmp[r,], show_rownames = F, show_colnames = F, annotation_col = geo_ok[,c("V11", "V7")])

head(pca2$u[,1:2])

plot(sample_af_means)
delta <- (sample_af_means[,1] - sample_af_means[,2])
plot(delta)

good <- names(tail(sort(abs(delta))))
par(mfrow = c(3,2))
lapply(good, function(x) boxplot(list(af[x, sample(asplit[[1]], 100)], af[x, sample(asplit[[2]], 100)]), main = x))

# Define (non-)infiltrating cells ####
subm <- meta[which(meta$Cluster_2d %in% c(11, 1)),]
subm2 <- subm[which(subm$Sample.name == "BT_S2"),]
groups <- split(rownames(subm2), subm2$Location)
subm <- meta[which(meta$Cluster_2d %in% c(11, 1)),]
infiltrating_one <- rownames(subm)[which(subm$Location == "Periphery")]
not_infiltrating_one <- rownames(subm)[-which(subm$Location == "Periphery")]

geo <- read.delim("/Users/lukas.simon/OneDrive/Miko/UTHealth/MitoTrace/DarmanisEtAl/FromGeoSpreadsheet.tsv", row.names = 1)
subm <- geo[which(geo$V11 == "Neoplastic"),]
subm2 <- subm[which(subm$V8 == "BT_S2"),]
infiltrating_one <- rownames(subm2)[which(subm2$V7 == "Tumor")]
not_infiltrating_one <- rownames(subm2)[which(subm2$V7 == "Periphery")]

# Load cell info ####
info <- read.delim("/Users/lukas.simon/OneDrive/Miko/UTHealth/MitoTrace/DarmanisEtAl/GBM_SRA_plate_relation")
cells <- c(infiltrating_one, not_infiltrating_one)
treat <- c(rep("inf", length(infiltrating_one)), rep("not", length(not_infiltrating_one)))

# Subset data ####
IDs <- as.character(info$sra[match(cells, info$plate.well)])
af_ok <- data.matrix(af[,IDs])
cov_ok <- data.matrix(cov[,IDs])
af_ok[which(cov_ok < 5)] <- NA
rownames(cov_ok) <- rownames(af_ok) <- rownames(af)

# Restrict to SNPs with min allele freq > 0.1 in either inf or not ####
asplit <- split(IDs, treat)
means <- do.call(cbind, lapply(asplit, function(x) rowMeans(af_ok[,x], na.rm = T)))
ok <- names(which(apply(means, 1, function(x) sum(x > 0.1)) > 0))

plot(means[,1], means[,2], xlab = "Infiltrating cells", ylab = "Non-infiltrating cells")
which(means[,1] > 0.2 & means[,2] == 0)

# Calculate differential allele freq between inf and not ####
pvals <- lapply(ok, function(x){
  print(x)
  af_tmp <- af_ok[x,]
  cov_tmp <- cov_ok[which(rownames(af_ok) == x),]
  af_tmp[which(cov_tmp <= 10)] <- NA
  group_inf <- split(af_tmp, treat)[["inf"]]
  group_not <- split(af_tmp, treat)[["not"]]
  try(wilcox.test(group_not, group_inf, na.rm=T)$p.value)
})
names(pvals) <- ok
bad <- which(unlist(lapply(pvals, class)) == "try-error")
if(length(bad) > 0) pvals <- pvals[-bad]
sort(unlist(pvals))

# Generate heatmap ####
tmp <- af_ok[ok,]
anno_col <- data.frame(treat, patient = meta[cells, "Sample.name"])
rownames(anno_col) <- colnames(tmp)
pheatmap::pheatmap(tmp, cluster_cols = T, annotation_col = anno_col, show_colnames = F)
