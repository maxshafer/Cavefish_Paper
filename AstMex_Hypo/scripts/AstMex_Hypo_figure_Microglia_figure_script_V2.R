library(Seurat)
library(stringr)
library(dplyr)
library(data.table)
library(purrr)
library(patchwork)
library(viridis)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")
hypo.ast <- readRDS("AstMex_63k_vR.rds")

# Subset out the blood lineage cells
Idents(hypo.ast) <- "Cluster"
immune <- subset(hypo.ast, idents = c("Erythrocytes", "Tcells", "Bcells", "Mast_cells", "Neutrophils", "Macrophages", "Microglia"))
immune <- FindVariableFeatures(immune, selection.method = "mvp")
immune <- ScaleData(object = immune, features = VariableFeatures(immune))
immune <- RunPCA(object = immune, features = VariableFeatures(immune), npcs = 100, set.seed = 0)
# ElbowPlot(object = immune, ndims = 100) # 25 PCs looks good
immune <- RunTSNE(object = immune, reduction = "pca", dims = 1:25, tsne.method = "Rtsne", reduction.name = "tsne", reduction.key = "tsne_", seed.use = 1, check_duplicates = F)

saveRDS(immune, file = "AstMex_immune.rds")

## Find Markers

cell.types <- unique(immune@meta.data$morph_Cluster)
cell.types <- cell.types[table(Idents(immune)) > 3]
Idents(immune) <- "morph_Cluster"
# Idents(immune) <- factor(Idents(immune), levels = levels(immune@meta.data$morph_Cluster))
morph_subtype_markers <- FindAllMarkers(immune, max.cells.per.ident = 500, only.pos = T)
markers <- morph_subtype_markers %>% group_by(cluster) %>% top_n(3, avg_logFC)


######### FILTER BY TRINARIZED GENES BEFORE PLOTTING #############

## Make Plots

cols0 <- c("#FDE725FF", "#22A884FF")
cols3 <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


immune.morph <- DimPlot(object = immune, group.by = "species", reduction = "tsne", pt.size = .25, label = FALSE, cols = cols0) + NoAxes() + theme(legend.position = c(0.8,0.9), legend.background = element_blank()) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
immune.subtype <- DimPlot(object = immune, group.by = "Cluster", reduction = "tsne", pt.size = .25, label = TRUE) + NoLegend() + NoAxes()
immune.subcluster <- DimPlot(object = immune, group.by = "Subcluster", reduction = "tsne", pt.size = .25, label = TRUE) + NoLegend() + NoAxes()
immune.orig <- DimPlot(object = immune, group.by = "orig.ident", reduction = "tsne", pt.size = .25, label = FALSE, cols = cols3) + NoAxes() + theme(legend.position = c(0.8,0.9), legend.background = element_blank()) + guides(color = guide_legend(ncol = 2, override.aes = list(size = 5))) + scale_colour_manual(values = cols3)

ccr9a <- FeaturePlot(immune, features = c("ccr9a"), reduction = "tsne", pt.size = .25) + NoAxes() + ggtitle("")

dot.plot <- DotPlot(immune, features = unique(markers$gene), group.by = "species_Cluster", scale.max = 150) + coord_flip() + theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_color_viridis(option = "A")

# Patchwork them together

tsnes <- (immune.subtype + immune.morph + ccr9a + plot_layout(nrow = 3)) 

tsnes / dot.plot + plot_layout(ncol = 2, widths = c(1.5,1), guides = "collect")



# ## Other gene lists, including from Peuss et. al.
# 
# gen.genes <- c("pfn1", "cd74a", "npc2", "grn2", "cotl1", "pgd")
# myeloid.genes <- c("apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b")
# lymphoid.genes <- c("cd28l", "sla2", "IGKC", "srgn", "rel", "p2ry11", "bric2", "ltb4r", "alox5ap")
# 
# 
# genes <- c("pfn1", "cd74a", "npc2", "grn2", "cotl1", "pgd", "apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b", "cd28l", "sla2", "IGKC", "srgn", "rel", "p2ry11", "bric2", "ltb4r", "alox5ap")
# 
# peuss.surface <- read.csv("~/Downloads/media-3_surface/Surface-Table 1.csv")
# peuss.pachon <- read.csv("~/Downloads/media-3_surface/Pachón-Table 1.csv")
# 
# top.pachon <- peuss.pachon %>% group_by(Cluster) %>% top_n(2, Avgerage.logFC)
# top.surface <- peuss.surface %>% group_by(Cluster) %>% top_n(2, Avgerage.logFC)
# pachon.genes <- top.pachon$Gene[top.pachon$Gene %in% row.names(GetAssayData(immune))]
# surface.genes <- top.surface$Gene[top.surface$Gene %in% row.names(GetAssayData(immune))]
# 
# DotPlot(object = immune, features = union(surface.genes, pachon.genes), group.by = "Subcluster") + RotatedAxis() + NoLegend()
# DotPlot(object = immune, features = genes, group.by = "Cluster") + RotatedAxis() + NoLegend()
# 
# print(top.pachon, n = 28)
# print(top.surface, n = 18)