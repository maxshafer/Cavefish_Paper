

## Need to set gene sets for heatmap: 1) Marker genes for each cluster, 2) Cave-specific and surface specific signatures, 3) Most Diff expressed genes, all cell types

## Plot the two tsne's (cave origin, plus cell type) - not including Erythrocytes

## Plot heatmap with genes signatures, like Bushra did, to show the cave signature, plus the individual highly DE genes


library(Seurat)
library(data.table)
library(stringr)
library(dplyr)
library(data.table)
library(purrr)
library(ggrepel)
library(biomaRt)
library(patchwork)
library(SuperExactTest)

load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_65k.Robj")


## Fix idents (change leucocyte Idents)

gen.genes <- c("pfn1", "cd74a", "npc2", "grn2", "cotl1", "pgd")
myeloid.genes <- c("apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b")
lymphoid.genes <- c("cd28l", "sla2", "IGKC", "srgn", "rel", "p2ry11", "bric2", "ltb4r", "alox5ap")


genes <- c("pfn1", "cd74a", "npc2", "grn2", "cotl1", "pgd", "apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b", "cd28l", "sla2", "IGKC", "srgn", "rel", "p2ry11", "bric2", "ltb4r", "alox5ap")



peuss.surface <- read.csv("~/Downloads/media-3_surface/Surface-Table 1.csv")
peuss.pachon <- read.csv("~/Downloads/media-3_surface/Pachón-Table 1.csv")

top.pachon <- peuss.pachon %>% group_by(Cluster) %>% top_n(2, Avgerage.logFC)
top.surface <- peuss.surface %>% group_by(Cluster) %>% top_n(2, Avgerage.logFC)
pachon.genes <- top.pachon$Gene[top.pachon$Gene %in% row.names(immune@data)]
surface.genes <- top.surface$Gene[top.surface$Gene %in% row.names(immune@data)]

dev.new()
DotPlot(object = immune, genes.plot = union(surface.genes, pachon.genes), plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType")

print(top.pachon, n = 28)
print(top.surface, n = 18)




immune <- SubsetData(hypo.ast, ident.use = c("Immune", "Blood"))
immune <- FindVariableGenes(immune, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE)
immune <- ScaleData(immune, genes.use = immune@var.genes, do.par = TRUE, num.cores = 6)
immune <- RunPCA(immune, pc.genes = immune@var.genes, pcs.compute = 25)
immune <- RunTSNE(immune, reduction.use = "pca", dims.use = 1:15, nthreads = 6, do.fast = TRUE, check_duplicates = FALSE)

cols3 <- c("goldenrod1", "goldenrod1", "lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1", "darkorange1", "darkorange1", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4"," springgreen4") # Individal caves/surface

DimPlot(object = immune, group.by = "orig.ident", reduction.use = "tsne", pt.size = .5, do.label = FALSE, no.legend = F, no.axes = T, cols.use = cols3)
DimPlot(object = immune, group.by = "species", reduction.use = "tsne", pt.size = .5, do.label = FALSE, no.legend = T, no.axes = T, cols.use = c("lightgoldenrod1", "springgreen4"))
DimPlot(object = immune, group.by = "SubclusterType", reduction.use = "tsne", pt.size = .5, do.label = TRUE, no.legend = TRUE, no.axes = T)
DimPlot(object = immune, group.by = "Subtype", reduction.use = "tsne", pt.size = .5, do.label = TRUE, no.legend = TRUE, no.axes = T)



test <- FeaturePlot(immune, features.plot = c("ccr9a"), cols.use = c("grey85", "blue"), no.axes = T, do.return = T, no.legend = F) + theme(text = element_blank())

## Find morph_subtype marker genes for plotting (maybe works best)

immune@meta.data$morph_subtype <- paste(immune@meta.data$species_morph, immune@meta.data$Subtype, sep = "_")

immune@meta.data$subtype_species <- paste(immune@meta.data$Subtype, immune@meta.data$species, sep = "_")

immune@meta.data$subtype_morph <- paste(immune@meta.data$Subtype, immune@meta.data$species_morph, sep = "_")


cell.types <- c("Bcells", "Mast_cells", "Tcells", "Microglia", "Macrophages", "Erythrocytes", "Thrombocytes", "Neutrophils")
cell.types <- unique(hypo.ast@meta.data$morph_subtype)[unlist(lapply(cell.types, function(x) grep(x, unique(hypo.ast@meta.data$morph_subtype))))]

hypo.ast <- SetAllIdent(hypo.ast, id = "morph_subtype")

morph_subtype_markers <- list()
morph_subtype_markers <- lapply(cell.types, function(x) FindMarkers(immune, ident.1 = x, max.cells.per.ident = 500, only.pos= T))

morph_subtype_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

markers$cluster <- immune@meta.data$subtype_morph[match(markers$cluster, immune@meta.data$morph_subtype)]

DotPlot(immune, genes.plot = unique(markers$gene), group.by = "subtype_species", do.return = T) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8)) + scale_color_viridis() + coord_flip()

plots <- lapply(cell.types, function(x) DotPlot(immune, genes.plot = unique(markers$gene[grep(x, markers$cluster)]), group.by = "subtype_species", do.return = T) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8)))
