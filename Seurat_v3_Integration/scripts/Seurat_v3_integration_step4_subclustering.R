library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")
load("Hypo_integrated_130k_1500VFs_100Dims_v0.Robj")
subsets <- readRDS("Hypo_integrated_130k_1500VFs_100Dims_subsets.rds")
neuronal <- subsets[[1]]

DefaultAssay(neuronal) <- "integrated"
neuronal <- FindNeighbors(neuronal)
neuronal <- FindClusters(neuronal, reduction.type = "pca", dims.use = 1:100, resolution = 0.2, print.output = 0)
DimPlot(neuronal, reduction = "tsne", group.by = "integrated_snn_res.0.2", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()

Idents(neuronal) <- "integrated_snn_res.0.35"
DefaultAssay(neuronal) <- "RNA"

DotPlot(neuronal, features = c("gng3", "slc17a6a", "gad1b", "slc32a1"), group.by = "integrated_snn_res.0.35") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

neuronal <- RenameIdents(neuronal, '0' = "GABA_0", '1' = "GABA_1", '2' = "Glut_0", '3' = "Glut_1", '4' = "Glut_2", '5' = "GABA_2", '6' = "Glut_3", '7' = "Glut_4", '8' = "Glut_5", '9' = "GABA_2", '10' = "Glut_6", '11' = "GABA_3", '12' = "GABA_4", '13' = "GABA_5", '14' = "GABA_6")

DimPlot(neuronal, reduction = "tsne", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()


# int.idents plus specific number of dimensions for integration/PCA
int.idents <- levels(Idents(neuronal))
int.dims <- c(50, 25, 20, 20, 20, 20, 15, 15, 10, 10, 5, 5, 5, 5)

# For each integrated subtype, re-run CCA/MNN integration
subsets.neuronal <- list()
for (i in 1:length(int.idents)) {
	subset <- subset(neuronal, idents = int.idents[[i]])
	Idents(subset) <- "species.2"
	reference.list <- c(subset(subset, idents = "zebrafish"), subset(subset, idents = "astyanax"))
	for (j in 1:length(reference.list)) {
		reference.list[[j]] <- NormalizeData(object = reference.list[[j]], verbose = FALSE)
		reference.list[[j]] <- FindVariableFeatures(object = reference.list[[j]], selection.method = "vst", nfeatures = 1500, verbose = FALSE)
	}
	# Integrated subsets
	subset.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:int.dims[[i]])
	subset.integrated <- IntegrateData(anchorset = subset.anchors, dims = 1:int.dims[[i]])
	# Scale and run dimensionality reduction on integrated data
	DefaultAssay(object = subset.integrated) <- "integrated"
	subset.integrated <- ScaleData(object = subset.integrated, verbose = FALSE)
	subset.integrated <- RunPCA(object = subset.integrated, npcs = int.dims[[i]], verbose = FALSE)
	subset.integrated <- RunTSNE(object = subset.integrated, reduction = "pca", dims = 1:int.dims[[i]], check_duplicates = F)
	subsets.neuronal[[i]] <- subset.integrated
}

names(subsets.neuronal) <- int.idents

# Combine into one object (neuronal and non-neuronal clusters - subclustered)
# Find subclusters

subsets.all <- c(subsets.neuronal, subsets[2:length(subsets)])

for (i in 1:length(subsets.all)) {
	DefaultAssay(subsets.all[[i]]) <- "integrated"
}

subsets.all <- lapply(subsets.all, function(x) FindNeighbors(x, reduction = "pca", dims = 1:length(x@reductions$pca)))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, reduction.type = "pca", dims.use = 1:length(x@reductions$pca), resolution = 0.25, print.output = 0))

# subsets.all[[10]]@meta.data$integrated_snn_res.0.4 <- 0
# subsets.all[[13]]@meta.data$integrated_snn_res.0.4 <- 0

# Use subsets.all to add meta.data to original hypo.integrated

ident.names <- names(subsets.all)

for (i in 1:length(subsets.all)) {
	hypo.integrated@meta.data[ WhichCells(subsets.all[[i]]) , "integrated_Subtype" ] <- ident.names[i]
}

# add subclustering labels
# for the rows of hypo.integrated meta data that match the cells of each subset (ident = integrated_snn_res.0.75), add $integrated_SubclusterType with paste(name of subset, ident number, sep = "_")

for (i in 1:length(subsets.all)) {
	Idents(subsets.all[[i]]) <- "integrated_snn_res.0.2"
}

for (i in 1:length(subsets.all)) {
	for (j in 1:length(levels(Idents(subsets.all[[i]])))) {
		hypo.integrated@meta.data[ WhichCells(subsets.all[[i]], idents = levels(Idents(subsets.all[[i]]))[j]) , "integrated_SubclusterType" ] <- paste(ident.names[i], levels(Idents(subsets.all[[i]]))[j], sep = "_")
	}
}

for (i in 1:length(subsets.all)) {
	for (j in 1:length(levels(Idents(subsets.all[[i]])))) {
		hypo.integrated@meta.data[ WhichCells(subsets.all[[i]], idents = levels(Idents(subsets.all[[i]]))[j]) , "integrated_SubclusterType_number" ] <- levels(Idents(subsets.all[[i]]))[j]
	}
}


for (i in 1:length(subsets.all)) {
	subsets.all[[i]]@meta.data$integrated_SubclusterType <- paste(ident.names[[i]], subsets.all[[i]]@meta.data$integrated_snn_res.0.2)
}

Idents(hypo.integrated) <- "integrated_Subtype"
levels(Idents(hypo.integrated))

# Remove GABA_5, which is a weird cluster of begotten cells, and rename remaining GABA clusters (do this for subclusters too, to save time)

hypo.integrated <- subset(hypo.integrated, cells = WhichCells(hypo.integrated, idents = c("GABA_5"), invert = T))

# rename

hypo.integrated <- RenameIdents(hypo.integrated, "GABA_6" = "GABA_5")

# Save files

save(hypo.integrated, file = "Hypo_integrated_130k_1500VFs_100Dims_v1.Robj")
saveRDS(subsets.all, file = "Hypo_integrated_130k_1500VFs_100Dims_subsets.rds")


## Remove weird subclusters (from species-specific analysis), but don't rename others
## Save as V2 (change old V2 as V1, and old V1 as V0)

Idents(hypo.integrated) <- "integrated_SubclusterType"

hypo.integrated <- subset(hypo.integrated, cells = WhichCells(hypo.integrated, idents = c("Glut_1_4", "Glut_3_7", "Glut_4_6", "GABA_2_9"), invert = T))
hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_4_7"] <- "Glut_4_6"

hypo.integrated@meta.data$integrated_SubclusterType <- factor(hypo.integrated@meta.data$integrated_SubclusterType, levels(hypo.integrated@meta.data$integrated_SubclusterType)[!(levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% c("Glut_1_4", "Glut_3_7", "Glut_4_7", "GABA_2_9"))])

hypo.integrated@meta.data$integrated_Subtype <- factor(hypo.integrated@meta.data$integrated_Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

save(hypo.integrated, file = "Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

