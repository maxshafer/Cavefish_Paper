library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v1.Robj")

hypo.integrated <- FindNeighbors(object = hypo.integrated)
hypo.integrated <- FindClusters(object = hypo.integrated, reduction.type = "pca", dims.use = 1:100, resolution = 0.2, print.output = 0)
hypo.integrated <- FindClusters(object = hypo.integrated, reduction.type = "pca", dims.use = 1:100, resolution = 0.35, print.output = 0)
hypo.integrated <- FindClusters(object = hypo.integrated, reduction.type = "pca", dims.use = 1:100, resolution = 0.8, print.output = 0)


cols <- c("#FFBE00", "#FF4500","#3E94D1") #Red/Orange/Blue

species <- DimPlot(object = hypo.integrated, reduction = "tsne", group.by = "species", pt.size = .1, cols = cols) + NoLegend() + NoAxes()
subtype <- DimPlot(object = hypo.integrated, reduction = "tsne", group.by = "Subtype", label = TRUE, repel = TRUE, pt.size = .5) + NoLegend() + NoAxes()

res.0.2 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.2", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()
res.0.35 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.35", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()
res.0.8 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.8", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()

plot_grid(subtype, res.0.2, res.0.35, res.0.8, nrow = 2)

DefaultAssay(hypo.integrated) <- "RNA"
FeaturePlot(hypo.integrated, features = c("cd74a"), reduction = "tsne", min.cutoff = "q9")


DotPlot(hypo.integrated, features = c("gng3", "slc17a6a", "slc1a2b", "gad1b", "slc6a1b", "slc32a1", "prdx1", "plp1b", "mpz", "timp4.1", "her15.1", "sfrp5", "vim", "arl13b", "pde6d", "hbaa2", "mb", "abcb4", "epd", "fetub", "mrc1a", "stab1", "pfn1", "cd74a", "npc2", "pgd", "cotl1", "tspan18a", "cxcr4b", "IGKC", "srgn", "rel", "ltb4r", "alox5ap", "dusp2", "apoeb", "apoc1", "cd28l", "sla2"), group.by = "integrated_snn_res.0.2") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

DotPlot(hypo.integrated, features = c("gng3", "gad1b", "galn", "oxt", "cd74b"), group.by = "integrated_snn_res.0.2")

Idents(hypo.integrated) <- "integrated_snn_res.0.35"

hypo.integrated <- RenameIdents(hypo.integrated, '0' = "Neuronal", '1' = "Neuronal", '2' = "Prdx1_Positive", '3' = "Neuronal", '4' = "Progenitors", '5' = "Leucocytes", '6' = "Neuronal", '7' = "Microglia", '8' = "Neuronal", '9' = "Oligodendrocytes", '10' = "Macrophages", '11' = "Endothelial", '12' = "Progenitors", '13' = "Oligodendrocyte_Precursor_Cells", '14' = "Ciliated", '15' = "Erythrocytes", '16' = "Lymphatic", '17' = "Ependymal", '18' = "Progenitors")

idents <- DimPlot(hypo.integrated, reduction = "tsne", pt.size = .1, label = T) + NoLegend() + NoAxes()

hypo.integrated$int.Subtype <- Idents(hypo.integrated)

plot_grid(idents, res.0.2, res.0.35, res.0.8, nrow = 2)


hypo.integrated@meta.data$species.2 <- "astyanax"
hypo.integrated@meta.data$species.2[hypo.integrated@meta.data$species == "zebrafish"] <- "zebrafish"

save(hypo.integrated, file = "Hypo_integrated_130k_1500VFs_100Dims_v1.Robj")

