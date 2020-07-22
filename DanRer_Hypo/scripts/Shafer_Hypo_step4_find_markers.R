library(Seurat)
library(Matrix)
library(dplyr)
library(data.table)
library(purrr)
library(data.table)
library(reshape2)


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo")


# Load marker genes
dir.markers <- paste("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/CSV/full_dataset/cluster_markers/", list.files("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/CSV/full_dataset/cluster_markers/"), sep = "")


markers.Subtype <- ldply(dir.markers, read.csv)

# Replace cluster numbers with names?
mappingvalues <- read.csv("CSV/full_dataset/Cluster_annotations_Shafer_Hypo.res.0.6.csv", header = TRUE, row.names = "cluster")

hypo <- SetAllIdent(hypo, id = "res.0.6")

current.cluster.ids <- c(0:43)

new.cluster.ids <- as.vector(mappingvalues$Subtype)

markers.Subtype$cluster <- mapvalues(x = markers.Subtype$cluster, from = current.cluster.ids, to = new.cluster.ids)

markers.Subtype <- markers.Subtype[order(match(markers.Subtype$cluster, levels(hypo@ident))),]

# Save Subtype markers

saveRDS(markers.Subtype, file = "Shafer_Hypo_markers.Subtype.rds")

markers.Subtype <- readRDS(file = "Shafer_Hypo_markers.Subtype.rds")

# Make figures

hypo <- SetAllIdent(hypo, id = "Subtype")

top2 <- markers.Subtype %>% group_by(cluster) %>% top_n(2, avg_logFC)
top10 <- markers.Subtype %>% group_by(cluster) %>% top_n(10, avg_logFC)
# trace(DotPlot, edit = t) to modify dot size, rotate axis before plotting below

## Load Seurat 3.0 to make proper DotPlots
detach("package:Seurat", unload = T)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
hypo <- UpdateSeuratObject(hypo.zeb)

rerio.marker.dot <- DotPlot(hypo, features = unique(top2$gene), group.by = "Subtype", dot.scale = 3, do.return = T) + coord_flip() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6), axis.title = element_blank()) + scale_color_viridis()

detach("package:Seurat", unload = T)
library(Seurat)



png("Figures/hypo_cluster_dotplot_Subtype_markers.png", height = 15, width = 8, units = "in", res = 250)
p2 <- DotPlot(object = hypo, genes.plot = unique(top2$gene), plot.legend = TRUE, do.return = T, dot.scale = 3) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6)) + scale_color_viridis() + coord_flip()
dev.off()

# Plot p2 with p1 from the AstMex file using patchwork
# p2 + p1 + plot_layout(guides = "collect", widths = c(length(unique(hypo.zeb@meta.data$Subtype)), length(unique(hypo.ast@meta.data$Subtype))))

# Legend on side again
png("Figures/hypo_cluster_dotplot_Subtype_markers_10.png", height = 40, width = 12, units = "in", res = 250)
DotPlot(object = hypo, genes.plot = unique(top10$gene[1:420]), plot.legend = TRUE, x.lab.rot = TRUE)
dev.off()


# Subcluster marker genes
dir.markers.sub <- paste("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/Figures/subclustering/CSV/", list.files("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/Figures/subclustering/CSV/"), sep = "")

subcluster.names <- list.files("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/Figures/subclustering/CSV/")
subcluster.names <- str_sub(subcluster.names, start = 14, end = -13)
subcluster.names <- mapvalues(x = subcluster.names, from = c(0:32), to = as.vector(mappingvalues $Subtype[1:33]))

markers.SubclusterType <- lapply(dir.markers.sub[c(1:26,28:33)], function(x) read.csv(x))
sub

names(markers.SubclusterType) <- subcluster.names[c(1:26,28:33)]


# save cluster markers
markers.SubclusterType <- melt(markers.SubclusterType, id.vars = c("cluster", "gene"), measure.vars = "p_val_adj")
markers.SubclusterType$SubclusterType <- paste(markers.SubclusterType$L1, markers.SubclusterType$cluster, sep = "_")

hypo <- SetAllIdent(hypo, id = "Subtype")

markers.SubclusterType <- markers.SubclusterType[order(match(markers.SubclusterType$L1, levels(hypo@ident))),]

# Save markers

saveRDS(markers.SubclusterType, file = "Shafer_Hypo_markers.SubclusterType.rds")

# Plot

top2 <- markers.SubclusterType %>% group_by(SubclusterType) %>% top_n(-2, value)

png("Figures/hypo_cluster_dotplot_SubclusterType_markers.png", height = 40, width = 20, units = "in", res = 250)
DotPlot(object = hypo, genes.plot = unique(top2$gene), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE)
dev.off()

