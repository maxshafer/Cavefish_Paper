library(Seurat)

## This is done on the server, scp the Seurat objects to the cluster, and run the following with SLURM
## Used R 4.0.0 on server

setwd("/scicore/home/schiera/gizevo30/projects/Seurat_v3_Integration")

hypo.integrated <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_v1.rds")

Idents(hypo.integrated) <- "integrated_Subtype"

# int.idents plus specific number of dimensions for integration/PCA
int.idents <- levels(Idents(hypo.integrated))
int.dims <- c(100, 50, 50, 35, 30, 25, 25, 20, 15, 15, 10, 10, 10) # dims1
# int.dims <- c(50, 25, 25, 15, 15, 10, 10, 10, 10, 7, 7, 7, 7) # dims2

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
	subsets[[i]] <- subset.integrated
}

# # Scale and run dimensionality reduction on integrated data
for (i in 1:length(subsets)) {
	DefaultAssay(subsets[[i]]) <- "integrated"
	subsets[[i]] <- ScaleData(subsets[[i]], verbose = T)
}

for (i in 1:length(subsets)) {
	subsets[[i]] <- RunPCA(subsets[[i]], npcs = int.dims[[i]], verbose = T)
}

for (i in 1:length(subsets)) {
	subsets[[i]] <- RunTSNE(subsets[[i]], reduction = "pca", dims = 1:int.dims[[i]], check_duplicates = F)
}

names(subsets) <- int.idents

saveRDS(subsets, file = "Hypo_integrated_128k_1500VFs_100Dims_subsets_dims1.rds")
