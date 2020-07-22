library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("DanRer_66k.Robj") # Loads as hypo

hypo.zeb <- UpdateSeuratObject(hypo.zeb)

load("AstMex_64k.Robj") # loads as hypo

hypo.ast <- UpdateSeuratObject(hypo)

rm(hypo)

reference.list <- c(hypo.zeb, hypo.ast)

for (i in 1:length(x = reference.list)) {
	Idents(reference.list[[i]]) <- reference.list[[i]]@meta.data$orig.ident
	reference.list[[i]] <- SubsetData(object = reference.list[[i]], max.cells.per.ident = 2000)
	# reference.list[[i]] <- SubsetData(object = reference.list[[i]], subset.name = "neuronal", accept.value = "neuronal")
    reference.list[[i]] <- NormalizeData(object = reference.list[[i]], verbose = FALSE)
    reference.list[[i]] <- FindVariableFeatures(object = reference.list[[i]], 
        selection.method = "vst", nfeatures = 1500, verbose = FALSE)
}

hypo.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:50)

hypo.integrated <- IntegrateData(anchorset = hypo.anchors, dims = 1:50)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData

DefaultAssay(object = hypo.integrated) <- "integrated"


# Run the standard workflow for visualization and clustering

hypo.integrated <- ScaleData(object = hypo.integrated, verbose = FALSE)
hypo.integrated <- RunPCA(object = hypo.integrated, npcs = 100, verbose = FALSE)
# hypo.integrated <- RunUMAP(object = hypo.integrated, reduction = "pca", dims = 1:100)
hypo.integrated <- RunTSNE(object = hypo.integrated, reduction = "pca", dims = 1:100, check_duplicates = F)


hypo.integrated.2 <- NormalizeData(hypo.integrated, block.size = 25000)
hypo.integrated.2 <- FindVariableFeatures(hypo.integrated.2, nfeatures = 2000)
hypo.integrated.2 <- ScaleData(object = hypo.integrated.2, features = VariableFeatures(hypo.integrated.2), verbose = FALSE)
hypo.integrated.2 <- RunPCA(object = hypo.integrated.2, assay = "RNA", reduction.key = "base", reduction.name = "base", npcs = 100, verbose = FALSE)
hypo.integrated.2 <- RunTSNE(object = hypo.integrated.2, reduction = "base", reduction.key = "basetsne_", reduction.name = "basetsne", dims = 1:100, check_duplicates = F)


save(hypo.integrated, file = "Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")