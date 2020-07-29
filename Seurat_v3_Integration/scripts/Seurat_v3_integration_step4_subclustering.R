library(Seurat)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_v1.rds")

subsets <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_subsets.rds")

neuronal <- subsets[[1]]

DefaultAssay(neuronal) <- "integrated"
neuronal <- FindNeighbors(neuronal, dims = 1:50, k.param = 20, prune.SNN = 1/15, nn.eps = 0.5)
neuronal <- FindClusters(neuronal, resolution = 0.1) # 15 clusters ~ similar to first clustering
neuronal <- FindClusters(neuronal, resolution = 0.2)
neuronal <- FindClusters(neuronal, resolution = 0.3)
DimPlot(neuronal, reduction = "tsne", group.by = "integrated_snn_res.0.1", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()

Idents(neuronal) <- "integrated_snn_res.0.1"
DefaultAssay(neuronal) <- "RNA"

DotPlot(neuronal, features = c("gng3", "slc17a6a", "gad1b", "slc32a1")) + RotatedAxis()

neuronal <- RenameIdents(neuronal, '0' = "GABA_0", '1' = "Glut_0", '2' = "Glut_1", '3' = "Glut_2", '4' = "Glut_3",
                                    '5' = "GABA_1", '6' = "Glut_4", '7' = "GABA_2", '8' = "GABA_3", '9' = "Glut_5", 
                                    '10' = "GABA_4", '11' = "GABA_5", '12' = "GABA_6", '13' = "Glut_6", '14' = "GABA_7")

DimPlot(neuronal, reduction = "tsne", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()


# int.idents plus specific number of dimensions for integration/PCA
int.idents <- levels(Idents(neuronal))
int.dims <- c(25, 15, 15, 15, 15, 15, 10, 10, 10, 10, 5, 5, 5, 5, 5)

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
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.25, print.output = 0))

# subsets.all[[10]]@meta.data$integrated_snn_res.0.4 <- 0
# subsets.all[[13]]@meta.data$integrated_snn_res.0.4 <- 0

# Use subsets.all to add meta.data to original hypo.integrated

ident.names <- names(subsets.all)

hypo.integrated@meta.data$integrated_Subtype <- factor(hypo.integrated@meta.data$integrated_Subtype, levels = c("Endothelial", "Erythrocytes", "Cilliated", "Ependymal", "Progenitors", "OPCs", "Oligodendrocytes", "Prdx1_Positive", sort(names(subsets.neuronal)), "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

for (i in 1:length(subsets.all)) {
	hypo.integrated@meta.data$integrated_Subtype[ row.names(hypo.integrated@meta.data) %in% row.names(subsets.all[[i]]@meta.data)] <- ident.names[[i]]
}

# add subclustering labels
# for the rows of hypo.integrated meta data that match the cells of each subset (ident = integrated_snn_res.0.75), add $integrated_SubclusterType with paste(name of subset, ident number, sep = "_")

for (i in 1:length(subsets.all)) {
	Idents(subsets.all[[i]]) <- "integrated_snn_res.0.25"
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
	subsets.all[[i]]@meta.data$integrated_SubclusterType <- paste(ident.names[[i]], subsets.all[[i]]@meta.data$integrated_snn_res.0.25)
}

Idents(hypo.integrated) <- "integrated_Subtype"
levels(Idents(hypo.integrated))

# Set factor levels for subclustertypes

index <- as.data.frame(unique(hypo.integrated@meta.data[,c("integrated_Subtype", "integrated_SubclusterType", "integrated_SubclusterType_number")]))
index <- index[order(index[,1], index[,3]),]

hypo.integrated@meta.data$integrated_SubclusterType <- factor(hypo.integrated@meta.data$integrated_SubclusterType, levels = index$integrated_SubclusterType)

# Save files

saveRDS(hypo.integrated, file = "Hypo_integrated_128k_1500VFs_100Dims_v1.rds")
saveRDS(subsets.all, file = "Hypo_integrated_128k_1500VFs_100Dims_subsets.rds")


## Remove weird subclusters (from species-specific analysis), but don't rename others
## Save as V2 (change old V2 as V1, and old V1 as V0)

Idents(hypo.integrated) <- "integrated_SubclusterType"

hypo.integrated <- subset(hypo.integrated, cells = WhichCells(hypo.integrated, idents = c("Glut_1_4", "Glut_3_7", "Glut_4_6", "GABA_2_9"), invert = T))
hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_4_7"] <- "Glut_4_6"

hypo.integrated@meta.data$integrated_SubclusterType <- factor(hypo.integrated@meta.data$integrated_SubclusterType, levels(hypo.integrated@meta.data$integrated_SubclusterType)[!(levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% c("Glut_1_4", "Glut_3_7", "Glut_4_7", "GABA_2_9"))])

hypo.integrated@meta.data$integrated_Subtype <- factor(hypo.integrated@meta.data$integrated_Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

save(hypo.integrated, file = "Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

