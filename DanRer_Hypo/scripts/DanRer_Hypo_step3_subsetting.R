library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

## Load objects

hypo <- readRDS("DanRer_65k.rds")

Idents(hypo) <- "Subtype"
idents <- levels(hypo@meta.data$Subtype)

# For each Subtype, re-run VariableFeatures, FindNeighbours, and FindClusters

# Subset by cluster, then re-cluster to identify subclusters
subsets <- lapply(idents, function(x) subset(hypo, idents = x))
names(subsets) <- idents

subsets <- lapply(subsets, function(x) NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000))
subsets <- lapply(subsets, function(x) FindVariableFeatures(x, selection.method = "mvp"))
subsets <- lapply(subsets, function(x) ScaleData(x, features = VariableFeatures(x)))
subsets <- lapply(subsets, function(x) RunPCA(object = x, features = VariableFeatures(x), npcs = 15, set.seed = 0))
subsets <- lapply(subsets, function(x) FindNeighbors(x, dims = 1:15, k.param = 30, nn.eps = 0.5))
subsets <- lapply(subsets, function(x) FindClusters(x, resolution = 0.4, random.seed = 0))
subsets <- lapply(subsets, function(x) FindClusters(x, resolution = 0.2, random.seed = 0))
subsets <- lapply(subsets, function(x) FindClusters(x, resolution = 0.6, random.seed = 0))



subsets <- lapply(subsets, function(x) RunTSNE(object = x, reduction = "pca", dims = 1:15, tsne.method = "Rtsne", reduction.name = "tsne", reduction.key = "tsne_", seed.use = 1, check_duplicates = F)) # Olig2 is too small for tsne

# Save subsets for further use
saveRDS(subsets, file = "DanRer_65k_subsets.rds")

# Add meta data to original object

for (i in 1:length(subsets)) {
  Idents(subsets[[i]]) <- "RNA_snn_res.0.4"
}

meta.data <- Reduce(rbind, lapply(subsets, function(x) x@meta.data))
meta.data$SubclusterType <- paste(meta.data$Subtype, meta.data$RNA_snn_res.0.4, sep = "_")
meta.data$SubclusterType_number <- meta.data$RNA_snn_res.0.4


hypo@meta.data$SubclusterType <- meta.data$SubclusterType[match(row.names(hypo@meta.data), row.names(meta.data))]
hypo@meta.data$SubclusterType_number <- meta.data$SubclusterType_number[match(row.names(hypo@meta.data), row.names(meta.data))]

Idents(hypo) <- "SubclusterType"

## Set factor levels for SubclusterTypes based on Subtype factors levels

index <- as.data.frame(unique(hypo@meta.data[,c("Subtype", "SubclusterType", "SubclusterType_number")]))
index <- index[order(index[,1], index[,3]),]

hypo@meta.data$SubclusterType <- factor(hypo@meta.data$SubclusterType, levels = index$SubclusterType)


saveRDS(hypo, file = "DanRer_65k.rds")
