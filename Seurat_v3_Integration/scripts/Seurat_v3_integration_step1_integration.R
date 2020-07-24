library(Seurat)

## This is done on the server, scp the two Seurat objects to the cluster, and run the following with SLURM
## Used R 4.0.0 on server

setwd("/scicore/home/schiera/gizevo30/projects/Seurat_v3_Integration")

hypo.zeb <- readRDS("DanRer_65k.rds")

hypo.ast <- readRDS("AstMex_63k.rds")

reference.list <- c(hypo.zeb, hypo.ast)

for (i in 1:length(x = reference.list)) {
  reference.list[[i]] <- NormalizeData(object = reference.list[[i]], verbose = FALSE)
  reference.list[[i]] <- FindVariableFeatures(object = reference.list[[i]], selection.method = "vst", nfeatures = 1500, verbose = FALSE)
}

hypo.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:100)

hypo.integrated <- IntegrateData(anchorset = hypo.anchors, dims = 1:100)

# switch to integrated assay. The variable features of this assay are
# automatically set during IntegrateData

DefaultAssay(object = hypo.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering

hypo.integrated <- ScaleData(object = hypo.integrated, verbose = FALSE)
hypo.integrated <- RunPCA(object = hypo.integrated, npcs = 100, verbose = FALSE)
hypo.integrated <- RunTSNE(object = hypo.integrated, reduction = "pca", dims = 1:100, check_duplicates = F)

save(hypo.integrated, file = "Hypo_integrated_128k_1500VFs_100Dims_v1.Robj")