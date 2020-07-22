library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(purrr)
library(data.table)
library(reshape2)


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo")


# Load marker genes
dir.markers <- paste("/Volumes/Maxwell/R_Projects/AstMex_Hypo/CSV/full_dataset/cluster_markers/", list.files("/Volumes/Maxwell/R_Projects/AstMex_Hypo/CSV/full_dataset/cluster_markers/"), sep = "")


markers.Subtype <- ldply(dir.markers, read.csv)

# Replace cluster numbers with names?
mappingvalues <- read.csv("CSV/full_dataset/Cluster_annotations_AstMex_Hypo.res.0.6.csv", header = TRUE, row.names = "cluster")

hypo <- SetAllIdent(hypo, id = "res.0.6")

current.cluster.ids <- as.numeric(levels(hypo@ident))

new.cluster.ids <- as.vector(mappingvalues$Subtype)[!grepl("CONT", mappingvalues$Subtype)]

markers.Subtype$cluster <- mapvalues(x = markers.Subtype$cluster, from = current.cluster.ids, to = new.cluster.ids)

markers.Subtype <- markers.Subtype[order(match(markers.Subtype$cluster, levels(hypo@ident))),]

# Save Subtype markers

saveRDS(markers.Subtype, file = "AstMex_Hypo_markers.Subtype.rds")

markers.Subtype <- readRDS(file = "AstMex_Hypo_markers.Subtype.rds")

markers.Subtype$cluster <- markers.Subtype %>% pull(cluster) %>% plyr::mapvalues(., c("Leucocytes_1", "Leucocytes_2", "Leucocytes_3", "Leucocytes_4"), c("Bcells", "Mast_cells", "Thrombocytes", "Neutrophils"))

markers.Subtype <- markers.Subtype[!grepl("32", markers.Subtype$cluster),]

## Change marker names for immune cell types

markers.Subtype$cluster <- factor(markers.Subtype$cluster, levels = levels(hypo@meta.data$Subtype))

# Make figures

hypo <- SetAllIdent(hypo, id = "Subtype")

top2 <- markers.Subtype %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- markers.Subtype %>% group_by(cluster) %>% top_n(10, avg_logFC)

top2 <- top2[order(top2$cluster),]
top10 <- top10[order(top10$cluster),]

## Load Seurat 3.0 to make proper DotPlots
detach("package:Seurat", unload = T)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
hypo <- UpdateSeuratObject(hypo.ast)

astyanax.marker.dot <- DotPlot(hypo, features = unique(top2$gene), group.by = "Subtype", dot.scale = 3, do.return = T) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6), axis.title = element_blank()) + scale_color_viridis()

library(patchwork)
rerio.marker.dot + astyanax.marker.dot + plot_layout(guides = "collect")

detach("package:Seurat", unload = T)
library(Seurat)


png("Figures/hypo_cluster_dotplot_Subtype_markers.png", height = 15, width = 8, units = "in", res = 250)
p1 <- DotPlot(object = hypo, genes.plot = unique(top2$gene), plot.legend = TRUE, do.return = T, dot.scale = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6)) + scale_color_viridis() + coord_flip()
dev.off()

# Plot p2 with p1 from the AstMex file using patchwork
# p2 + p1 + plot_layout(guides = "collect", widths = c(length(unique(hypo.zeb@meta.data$Subtype)), length(unique(hypo.ast@meta.data$Subtype))))

widths = c(length(unique(hypo.zeb@meta.data$Subtype)), length(unique(hypo.ast@meta.data$Subtype)))

# Legend on side again
png("Figures/hypo_cluster_dotplot_Subtype_markers_10.png", height = 40, width = 12, units = "in", res = 250)
DotPlot(object = hypo, genes.plot = unique(top10$gene), plot.legend = TRUE, x.lab.rot = TRUE, do.return = T) + scale_color_viridis()
dev.off()


# Subcluster marker genes
## Find Marker genes (the ones I have are for non-combined subclusters (~250))
hypo <- SetAllIdent(hypo, id = "Subtype")

hypo.small <- SubsetData(hypo, max.cells.per.ident = 1000)

hypo.small <- SetAllIdent(hypo.small, id = "SubclusterType")

markers.SubclusterType <- FindAllMarkers(object = hypo.small, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.5, max.cells.per.ident = 250)

# Save markers

saveRDS(markers.SubclusterType, file = "Shafer_Hypo_markers.SubclusterType.rds")

# Plot

top2 <- markers.SubclusterType %>% group_by(cluster) %>% top_n(2, avg_logFC)

hypo <- SetAllIdent(hypo, id = "Subtype")

png("Figures/hypo_cluster_dotplot_SubclusterType_markers.png", height = 80, width = 30, units = "in", res = 250)
DotPlot(object = hypo, genes.plot = unique(top10$gene), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE)
dev.off()

DotPlot(object = SubsetData(hypo, ident.use = c("Glut_2")), genes.plot = unique(top10$gene[grep("Glut_2", top10$cluster)]), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE)
