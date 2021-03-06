library(Seurat)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_v1.rds")

subsets <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_subsets_dims1.rds")

neuronal <- subsets[[1]]

DefaultAssay(neuronal) <- "integrated"
neuronal <- FindNeighbors(neuronal, dims = 1:10, k.param = 20, prune.SNN = 1/15, nn.eps = 0)
neuronal <- FindClusters(neuronal, resolution = 0.1) # 15 clusters ~ similar to first clustering
neuronal <- FindClusters(neuronal, resolution = 0.2)
neuronal <- FindClusters(neuronal, resolution = 0.3)
neuronal <- FindClusters(neuronal, resolution = 0.35)
DimPlot(neuronal, reduction = "tsne", group.by = "integrated_snn_res.0.35", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()

Idents(neuronal) <- "integrated_snn_res.0.35"
DefaultAssay(neuronal) <- "RNA"

DotPlot(neuronal, features = c("gng3", "slc17a6a", "slc17a6b", "gad2", "gad1b", "slc32a1")) + RotatedAxis()

neuronal <- RenameIdents(neuronal, '0' = "GABA_0", '1' = "GABA_1", '2' = "Glut_0", '3' = "Glut_1", '4' = "Glut_2",
                                    '5' = "GABA_2", '6' = "Glut_3", '7' = "Glut_4", '8' = "Glut_5", '9' = "GABA_3", 
                                    '10' = "GABA_4", '11' = "Glut_6", '12' = "Glut_7")

DimPlot(neuronal, reduction = "tsne", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()


# int.idents plus specific number of dimensions for integration/PCA
int.idents <- levels(Idents(neuronal))
# int.dims <- c(25, 15, 15, 15, 15, 15, 10, 10, 10, 10, 5, 5, 5, 5, 5)
int.dims <- c(50, 25, 25, 25, 25, 20, 20, 15, 15, 15, 10, 10, 5)

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
	subset.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:int.dims[[i]])  # add k.filter = 100 for int.idents[[15]]
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
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.2, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.25, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.4, print.output = 0))

#########################################################
## NEED TO ADD CILLIATED AS A SUBSET, and reorder
## REMOVE zeb cell? use classic subclustering (not integrated, but RNA assay)

subset <- subset(hypo.integrated, idents = "Cilliated")
# Idents(subset) <- "species.2"
# subset <- subset(subset, idents = "astyanax")
DefaultAssay(subset) <- "RNA"
subset <- NormalizeData(subset)
subset <- FindVariableFeatures(object = subset, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
subset <- ScaleData(subset, verbose = T)
subset <- RunPCA(subset, npcs = 10, verbose = T)
subset <- RunTSNE(subset, reduction = "pca", dims = 10, check_duplicates = F)
subset <- FindNeighbors(subset, reduction = "pca", dims = 1:length(subset@reductions$pca))
subset <- FindClusters(subset, resolution = 0.2, print.output = 0)
subset <- FindClusters(subset, resolution = 0.25, print.output = 0)
subset <- FindClusters(subset, resolution = 0.4, print.output = 0)
subset@meta.data$integrated_snn_res.0.25 <- subset@meta.data$RNA_snn_res.0.25
# Add back to subsets, and re-order
#########################################################

subsets.all.2 <- c(subsets.all, subset)
names(subsets.all.2) <- c(names(subsets.all), "Cilliated")
subsets.all <- subsets.all.2

# res 0.25 looks good, gives ~ 200 clusters

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
	subsets.all[[i]]@meta.data$integrated_SubclusterType <- paste(ident.names[[i]], subsets.all[[i]]@meta.data$integrated_snn_res.0.25, sep = "_")
}

Idents(hypo.integrated) <- "integrated_Subtype"
levels(Idents(hypo.integrated))

# Set factor levels for subclustertypes

index <- as.data.frame(unique(hypo.integrated@meta.data[,c("integrated_Subtype", "integrated_SubclusterType", "integrated_SubclusterType_number")]))
index <- index[order(index[,1], index[,3]),]
index <- unique(index[,c(1,2)])

hypo.integrated@meta.data$integrated_SubclusterType <- factor(hypo.integrated@meta.data$integrated_SubclusterType, levels = index$integrated_SubclusterType)

# Save files

saveRDS(hypo.integrated, file = "Hypo_integrated_128k_1500VFs_100Dims_v1.rds")
saveRDS(subsets.all, file = "Hypo_integrated_128k_1500VFs_100Dims_subsets_dims1.rds")

