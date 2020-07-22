library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration")

# Load whole integrated dataset and set idents to integrated subtypes
load("Hypo_integrated_130k_1500VFs_100Dims_v1.Robj")
Idents(hypo.integrated) <- "int.Subtype"

# int.idents plus specific number of dimensions for integration/PCA
int.idents <- levels(Idents(hypo.integrated))
int.dims <- c(100, 50, 50, 35, 30, 25, 25, 20, 15, 15, 10, 10, 10)

# For each integrated subtype, re-run CCA/MNN integration
subsets <- list()
for (i in 1:length(int.idents)) {
	subset <- subset(hypo.integrated, idents = int.idents[[i]])
	Idents(subset) <- "species.2"
	reference.list <- c(subset(subset, idents = "zebrafish"), subset(subset, idents = "astyanax"))
	for (j in 1:length(reference.list)) {
		reference.list[[j]] <- NormalizeData(object = reference.list[[j]], verbose = FALSE)
		reference.list[[j]] <- FindVariableFeatures(object = reference.list[[j]], selection.method = "vst", nfeatures = 1500, verbose = FALSE)
	}
	# Integrated subsets
	subset.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:int.dims[[i]])
	subset.integrated <- IntegrateData(anchorset = subset.anchors, dims = 1:int.dims[[i]], k.weight = min(100, min(ncol(reference.list[[1]]@assays$RNA), ncol(reference.list[[2]]@assays$RNA))))
	# Scale and run dimensionality reduction on integrated data
	DefaultAssay(object = subset.integrated) <- "integrated"
	subset.integrated <- ScaleData(object = subset.integrated, verbose = FALSE)
	subset.integrated <- RunPCA(object = subset.integrated, npcs = int.dims[[i]], verbose = FALSE)
	subset.integrated <- RunTSNE(object = subset.integrated, reduction = "pca", dims = 1:int.dims[[i]], check_duplicates = F)
	subsets[[i]] <- subset.integrated
}

names(subsets) <- int.idents

saveRDS(subsets, file = "Hypo_integrated_130k_1500VFs_100Dims_subsets.rds")
