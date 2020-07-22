###################################### Non-neuronal cells ######################################

hypo_neuronal <- CreateSeuratObject(input.data.nonneuronal, min.genes = 200, project = "Shafer_Hypo")
hypo_nonneuronal <- NormalizeData(object = hypo_nonneuronal, normalization.method = "LogNormalize", scale.factor = 10000)
hypo_nonneuronal <- ScaleData(object = hypo_nonneuronal)
hypo_nonneuronal <- FindVariableGenes(object = hypo_nonneuronal)
length(x = hypo_nonneuronal@var.genes)

hypo_nonneuronal <- RunPCA(object = hypo_nonneuronal, pc.genes = hypo_nonneuronal@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)

pdf("Figures/PCElbowPlot_hypo_nonneuronal.pdf")
PCElbowPlot(object = hypo_nonneuronal, num.pc = 100) # Looks like anywhere from 20-30
dev.off()

hypo_nonneuronal <- FindClusters(object = hypo_nonneuronal, reduction.type = "pca", dims.use = 1:25, resolution = 1.2, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

hypo_nonneuronal <- RunTSNE(object = hypo_nonneuronal, dims.use = 1:25, check_duplicates = FALSE, do.fast = TRUE)

pdf("Figures/hypo_nonneuronal_cluster_plot.pdf")
hypo_nonneuronal_plot <- TSNEPlot(object = hypo_nonneuronal, do.label = TRUE)
dev.off()

# FeaturePlot(object = hypo_nonneuronal, features.plot = c("sst1.2"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = .5)

# Find markers

hypo_nonneuronal.markers <- FindAllMarkers(object = hypo_nonneuronal, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# Generate Feature plot of top 2 markers each

features.markers <- hypo_nonneuronal.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

pdf("Figures/cluster.markers.hypo_nonneuronal.pdf")
FeaturePlot(object = hypo_nonneuronal, features.plot = features.markers$gene, min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf("Figures/cluster.markers.vlnplot.hypo_nonneuronal.pdf")
VlnPlot(object = hypo_nonneuronal, features.plot = features.markers$gene, point.size.use = NA)
dev.off()

# Generate lists + feature plots of top 20 markers for each cluster

for (i in 0:13) {
  name <- paste("cluster_nonneuronal_", i, sep="")
  assign(name, hypo_nonneuronal.markers %>% filter(cluster == i))
}

# Generate feature plots for each cluster (top 20 genes)

for (i in 0:13) {
  cluster <- hypo_nonneuronal.markers %>% filter(cluster == i) %>% top_n(10, avg_logFC)
  pdf(paste("Figures/hypo.cluster.markers.hypo_nonneuronal.", i, ".pdf", sep = ""))
  FeaturePlot(object = hypo_nonneuronal, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.1)
  dev.off()
}


top10 <- hypo_nonneuronal.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top20 <- hypo_nonneuronal.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

pdf("Figures/cluster.markers.hypo_nonneuronal.heatmap.top10.pdf")
DoHeatmap(hypo_nonneuronal, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
pdf("Figures/cluster.markers.hypo_nonneuronal.heatmap.top20.pdf")
DoHeatmap(hypo_nonneuronal, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()