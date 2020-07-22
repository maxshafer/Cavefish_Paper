

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

hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_1"] <- "Bcells"
hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_2"] <- "Mast_cells"
hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_3"] <- "Thrombocytes"
hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_4"] <- "Neutrophils"


# Load marker genes matrices and lists

# matrices
subtype.markers <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_markers.Subtype.rds")
subcluster.markers <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_markers.SubclusterType.rds")

# lists
subtype.markers.species <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_markers.species.Subtype.list.rds")
subcluster.markers.species <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_markers.species.SubclusterType.list.rds")

for (i in 1:length(subtype.markers.species)) {
	subtype.markers.species[[i]]$gene <- rownames(subtype.markers.species[[i]])
}

test <- lapply(subtype.markers.species, function(x) top_n(x, -5, p_val_adj))
test <- lapply(subtype.markers.species, function(x) rbind(top_n(x, 5, avg_logFC), top_n(x, -5, avg_logFC)))



## Extract gene lists

## SuperExactTest

zeb.markers.2 <- lapply(zeb.markers, function(x) x[x %in% zeb.genes])
ast.markers.2 <- lapply(ast.markers, function(x) x[x %in% ast.genes])

zeb.markers.3 <- zeb.markers[zeb.names]
ast.markers.3 <- ast.markers.2[ast.names]

length.gene.sets <- sapply(zeb.markers.3,length)

total <- nrow(hypo.zeb@data)

# Just neuronal zebrafish clusters
res=supertest(zeb.markers.3[c(2,3,5:9)], n=total)
	
plot(res, Layout="landscape", degree = 2:11, sort.by="size", keep=FALSE, show.elements=TRUE, elements.cex=0.7, show.fold.enrichment = TRUE, elements.list=subset(summary(res)$Table,Observed.Overlap <= 6), show.expected.overlap=TRUE,expected.overlap.style="hatchedBox", color.expected.overlap='red')





common.signature <- intersect(intersect(intersect(row.names(subtype.markers.species[["Macrophages"]]), row.names(subtype.markers.species[["Microglia"]])), intersect(row.names(subtype.markers.species[["Tcells"]]), row.names(subtype.markers.species[["Leucocytes_2"]]))), intersect(row.names(subtype.markers.species[["Leucocytes_3"]]), row.names(subtype.markers.species[["Leuococytes_4"]])))

common.signature <- intersect(intersect(row.names(subtype.markers.species[[26]]), row.names(subtype.markers.species[[27]])), intersect(row.names(subtype.markers.species[[36]]), row.names(subtype.markers.species[[21]])))

# Extract marker genes

top.subtype.markers <- subtype.markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
test <- top.subtype.markers[top.subtype.markers$cluster == c("Macrophages") | top.subtype.markers$cluster == c("Microglia") | top.subtype.markers$cluster == c("Tcells") | top.subtype.markers$cluster == c("Leucocytes_1") | top.subtype.markers$cluster == c("Leucocytes_2") | top.subtype.markers$cluster == c("Leucocytes_3") | top.subtype.markers$cluster == c("Leucocytes_4"),]

marker.genes <- as.character(test$gene)


("Microglia", "Tcells", "Leucocytes_1", "Leucocytes_2", "Leucocytes_3", "Leucocytes_4")


gen.genes <- c("pfn1", "cd74a", "npc2", "grn2", "cotl1", "pgd")
myeloid.genes <- c("apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b")
lymphoid.genes <- c("cd28l", "sla2", "IGKC", "srgn", "rel", "p2ry11", "bric2", "ltb4r", "alox5ap")


genes <- c("pfn1", "cd74a", "npc2", "grn2", "cotl1", "pgd", "apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b", "cd28l", "sla2", "IGKC", "srgn", "rel", "p2ry11", "bric2", "ltb4r", "alox5ap")

# Subset object for plotting heatmap

immune <- SubsetData(hypo.ast, ident.use = c("Immune"))
immune <- FindVariableGenes(immune, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE)
immune <- ScaleData(immune, do.par = TRUE, num.cores = 6)
immune <- RunPCA(immune, pc.genes = immune@var.genes, pcs.compute = 15)
immune <- RunTSNE(immune, reduction.use = "pca", dims.use = 1:15, nthreads = 6, do.fast = TRUE, check_duplicates = FALSE)

immune@meta.data$Subtype_species <- paste(immune@meta.data$Subtype, immune@meta.data$species, sep = "_")
immune@meta.data$Subtype_morph <- paste(immune@meta.data$Subtype, immune@meta.data$species_morph, sep = "_")


plot <- DoHeatmap(object = immune, cells.use = WhichCells(immune, max.cells.per.ident = 100), genes.use = c(gen.genes, marker.genes), group.by = "Subtype_species", group.order = NULL, title = NULL, slim.col.label = TRUE, group.label.rot = T)

plot <- DoHeatmap(object = immune, genes.use = c(test[[26]]$gene), group.by = "Subtype_morph", col.low = "#440154FF", col.mid = "#21908CFF", col.high = "#FDE725FF", group.order = NULL, title = NULL, slim.col.label = TRUE, group.label.rot = T)


plot <- DoHeatmap(immune, cells.use = WhichCells(immune, max.cells.per.ident = 20), genes.use = c(common.signature), group.by = "Subtype_species", group.order = NULL, title = NULL, slim.col.label = TRUE, group.label.rot = T)

# Plot with new color scale
plot + scale_fill_viridis(name = "Expression", guide = guide_colorbar(direction = key.direction, title.position = key.title.pos), option = "D")





### Load averaged data

norm.cluster <- readRDS("Normed_expression_data.rds")

str(norm.cluster, max.level = 2)

# Plot Surface vs cave log'd expression for microglia
plot(log(norm.cluster[[1]][[3]]$Progenitors$mean.exp), log(norm.cluster[[1]][[4]]$Progenitors$mean.exp))

norm.type2 <- do.call(cbind, c(norm.cluster[[1]][[3]], norm.cluster[[1]][[4]]))
colnames(norm.type2) <- c(paste("ast.surface", names(norm.cluster[[1]][[3]]), sep = "."), paste("ast.cave", names(norm.cluster[[1]][[4]]), sep = "."))
norm.type2 <- log(norm.type2)
norm.type2$gene <- row.names(norm.type2)

ScatterPlot <- function(data = norm.type2, cell.type = cell.type, subset.ids.1 = subset.ids, subset.ids.2 = subset.ids) {
	plot <- ggplot(data, aes_string(x = paste("ast.surface", cell.type, sep = "."), y = paste("ast.cave", cell.type, sep = ".")))
	plot <- plot + geom_point(size = 1, color = "grey85")
	
	plot <- plot + geom_text_repel(data = data[subset.ids.1,], aes_string(x = paste("ast.surface", cell.type, sep = "."), y = paste("ast.cave", cell.type, sep = "."), label = "gene"), color = "blue")
	plot <- plot + geom_point(data = data[subset.ids.1,], aes_string(x = paste("ast.surface", cell.type, sep = "."), y = paste("ast.cave", cell.type, sep = ".")), color = "blue")
	
	plot <- plot + geom_text_repel(data = data[subset.ids.2,], aes_string(x = paste("ast.surface", cell.type, sep = "."), y = paste("ast.cave", cell.type, sep = "."), label = "gene"), color = "blue")
	plot <- plot + geom_point(data = data[subset.ids.2,], aes_string(x = paste("ast.surface", cell.type, sep = "."), y = paste("ast.cave", cell.type, sep = ".")), color = "blue")
	
	return(plot)
}

ScatterPlot(data = norm.type2, cell.type = "Microglia", subset.ids.1 = common.signature, subset.ids.2 = common.signature)


#####
common.signature <- intersect(intersect(intersect(row.names(morph.markers[[26]]), row.names(morph.markers[[27]])), intersect(row.names(morph.markers[[36]]), row.names(morph.markers[[22]]))), intersect(row.names(morph.markers[[23]]), row.names(morph.markers[[24]])))

common.signature <- intersect(intersect(row.names(morph.markers[[26]]), row.names(morph.markers[[27]])), intersect(row.names(morph.markers[[36]]), row.names(morph.markers[[22]])))


DotPlot(object = SubsetData(hypo.ast, ident.use = "Immune"), genes.plot = common.signature, plot.legend = TRUE, x.lab.rot = TRUE, group.by = "species", do.return = T) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8))










immune <- SubsetData(hypo.ast, ident.use = c("Immune", "Blood"))
immune <- FindVariableGenes(immune, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE)
immune <- ScaleData(immune, genes.use = immune@var.genes, do.par = TRUE, num.cores = 6)
immune <- RunPCA(immune, pc.genes = immune@var.genes, pcs.compute = 25)
immune <- RunTSNE(immune, reduction.use = "pca", dims.use = 1:15, nthreads = 6, do.fast = TRUE, check_duplicates = FALSE)

cols3 <- c("goldenrod1", "goldenrod1", "lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1", "darkorange1", "darkorange1", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4"," springgreen4") # Individal caves/surface

DimPlot(object = immune, group.by = "orig.ident", reduction.use = "tsne", pt.size = .5, do.label = FALSE, no.legend = F, no.axes = T, cols.use = cols3)
DimPlot(object = immune, group.by = "species", reduction.use = "tsne", pt.size = .5, do.label = FALSE, no.legend = F, no.axes = T, cols.use = c("lightgoldenrod1", "springgreen4"))
DimPlot(object = immune, group.by = "SubclusterType", reduction.use = "tsne", pt.size = .5, do.label = TRUE, no.legend = TRUE, no.axes = T)
DimPlot(object = immune, group.by = "Subtype", reduction.use = "tsne", pt.size = .5, do.label = TRUE, no.legend = TRUE, no.axes = T)







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






DotPlot(object = immune, genes.plot = rownames(hypo.ast@data)[grep("sat", rownames(hypo.ast@data))], plot.legend = TRUE, x.lab.rot = TRUE, group.by = "species")



DotPlot(object = immune, genes.plot = astmex.markers[astmex.markers$cluster == "Leucocytes_2", "gene"][1:20], plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType")



## Find morph_subtype marker genes for plotting (maybe works best)

hypo.ast@meta.data$morph_subtype <- paste(hypo.ast@meta.data$species_morph, hypo.ast@meta.data$Subtype, sep = "_")

cell.types <- c("Bcells", "Mast_cells", "Tcells", "Microglia", "Macrophages", "Erythrocytes", "Thrombocytes", "Neutrophils")
cell.types <- unique(hypo.ast@meta.data$morph_subtype)[unlist(lapply(cell.types, function(x) grep(x, unique(hypo.ast@meta.data$morph_subtype))))]

hypo.ast <- SetAllIdent(hypo.ast, id = "morph_subtype")

morph_subtype_markers <- lapply(cell.types, function(x) FindMarkers(hypo.ast, ident.1 = x))



