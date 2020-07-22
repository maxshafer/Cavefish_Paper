library(Seurat)
library(Matrix)
library(dplyr)
library(datapasta)
library(corrgram)
library(corrplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo")

load("Shafer_Hypo_66k.Robj")


# Load GO lists

TF_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/GO_DNA_binding.csv", head = F)
TF_list <- as.character(unique(TF_list$V2))
TF_list <- TF_list[TF_list %in% rownames(hypo.zeb@data)]

NP_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/GO_neuropeptide.csv", head = F)
NP_list <- as.character(unique(NP_list$V2))
NP_list <- NP_list[NP_list %in% rownames(hypo.zeb@data)]


data <- as.matrix(hypo.zeb@data[NP_list,])
data[data > 1] <- 1
data[data < 1] <- 0

test <- colSums(data)
hist(test)

hypo.zeb@meta.data$peptide_count <- test


data <- as.matrix(hypo.zeb@data[TF_list,])
data[data < 2] <- 0
data[data > 2] <- 1

test <- colSums(data)
hist(test)

hypo.zeb@meta.data$tf_count <- test


DimPlot(object = hypo.zeb, group.by = "peptide_count", reduction.use = "tsne", pt.size = .05, do.label = F, label.size = 4, no.legend = F, no.axes = T)
DimPlot(object = hypo.zeb, group.by = "tf_count", reduction.use = "tsne", pt.size = .05, do.label = F, label.size = 4, no.legend = F, no.axes = T)
