library(Seurat)
library(ggplot2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_v1.rds")
# hypo.integrated <- readRDS("~/Documents/Seurat_objects/Hypo_integrated_128k_1500VFs_100Dims_v1.rds")

DefaultAssay(hypo.integrated) <- "integrated"

# Find clusters, and compare the cluster with the independent annotation

hypo.integrated <- FindNeighbors(object = hypo.integrated, dims = 1:100, k.param = 20, nn.eps = 0)
hypo.integrated <- FindClusters(object = hypo.integrated, resolution = 0.1)
hypo.integrated <- FindClusters(object = hypo.integrated, resolution = 0.15)
hypo.integrated <- FindClusters(object = hypo.integrated, resolution = 0.2)
hypo.integrated <- FindClusters(object = hypo.integrated, resolution = 0.3)

cols <- c("#FDE725FF", "#22A884FF", "#414487FF") # Viridis colours

species <- DimPlot(object = hypo.integrated, reduction = "tsne", group.by = "species", pt.size = .1, cols = cols) + NoLegend() + NoAxes()
subtype <- DimPlot(object = hypo.integrated, reduction = "tsne", group.by = "Subtype", label = TRUE, repel = TRUE, pt.size = .5) + NoLegend() + NoAxes()

res.0.1 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.1", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()
res.0.15 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.15", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()
res.0.2 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.2", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()
res.0.3 <- DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_snn_res.0.3", pt.size = .5, label = TRUE) + NoLegend() + NoAxes()

# It seems that at almost all resolutions, Tcells, Bcells and Mast_cells cluster together as a Leucocyte population
# At res.0.15, Macrophages and Microglia separate, and we get 11 cluster of neuronal cells
# We had best clustering results when we re-clustered the neuronal cells together, follow by subclustering
# We will therefore merge the neuronal clusters, and re-do the clustering on them alone, and recombine with the rest of the cells

subtype + res.0.1 + res.0.15 + res.0.2

# Check a few genes across clusters

DefaultAssay(hypo.integrated) <- "RNA"
Idents(hypo.integrated) <- "integrated_snn_res.0.15"

DotPlot(hypo.integrated, features = c("gng3", "slc17a6a", "slc1a2b", "gad1b", "slc6a1b", "slc32a1", "prdx1", "plp1b", "mpz", "timp4.1", "her15.1", "sfrp5", "vim", "arl13b", "pde6d", "hbaa2", "mb", "abcb4", "epd", "fetub", "mrc1a", "stab1", "pfn1", "cd74a", "npc2", "pgd", "cotl1", "tspan18a", "cxcr4b", "IGKC", "srgn", "rel", "ltb4r", "alox5ap", "dusp2", "apoeb", "apoc1", "cd28l", "sla2")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
DotPlot(hypo.integrated, features = c("gng3", "gad1b", "galn", "oxt", "cd74b"))

hypo.integrated <- RenameIdents(hypo.integrated, '0' = "Neuronal", '1' = "Neuronal", '2' = "Progenitors", '3' = "Prdx1_Positive", '4' = "Neuronal", 
                                                  '5' = "Leucocytes", '6' = "Neuronal", '7' = "Microglia", '8' = "Neuronal", '9' = "Macrophages", 
                                                  '10' = "Oligodendrocytes", '11' = "Neuronal", '12' = "Neuronal", '13' = "OPCs", 
                                                  '14' = "Endothelial", '15' = "Neuronal", '16' = "Cilliated", '17' = "Erythrocytes", '18' = "Lymphatic",
                                                  '19' = "Ependymal", '20' = "Neuronal", '21' = "Neuronal", '22' = "Neuronal", '23' = "Leucocytes")

hypo.integrated$integrated_Subtype <- Idents(hypo.integrated)


## Add meta data for zebrafish vs Mexican tetra
hypo.integrated@meta.data$species.2 <- "astyanax"
hypo.integrated@meta.data$species.2[hypo.integrated@meta.data$species == "zebrafish"] <- "zebrafish"

species.2 <- DimPlot(hypo.integrated, reduction = "tsne", pt.size = .1, group.by = "species.2") + NoLegend() + NoAxes()
idents <- DimPlot(hypo.integrated, reduction = "tsne", pt.size = .1, label = T) + NoLegend() + NoAxes()

species.2 + idents

saveRDS(hypo.integrated, file = "Hypo_integrated_128k_1500VFs_100Dims_v1.rds")

