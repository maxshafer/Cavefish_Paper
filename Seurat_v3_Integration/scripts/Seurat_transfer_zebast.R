library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/scicore/home/schiera/gizevo30/projects/Seurat_v3_Integration")

load("Hypo_integrated_130k_1500VFs_100Dims.Robj")

load("AstMex_66k.Robj")
hypo.ast <- UpdateSeuratObject(hypo)
rm(hypo)

Idents(hypo.integrated) <- "species"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")

hypo.anchors <- FindTransferAnchors(reference = hypo.integrated.zeb, query = hypo.ast, reduction = "pcaproject", dims = 1:100, npcs = NULL)

predictions <- TransferData(anchorset = hypo.anchors, refdata = hypo.zeb$Subtype, dims = 1:50)

hypo.ast <- AddMetaData(hypo.ast, metadata = predictions)

save(hypo.ast, file = "AstMex_66k_zeb.Robj")

