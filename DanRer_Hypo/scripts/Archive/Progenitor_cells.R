###################################### Progenitor cells ######################################

hypo_progenitors <- CreateSeuratObject(input.data.progenitors, min.genes = 200, project = "Shafer_Hypo")
hypo_progenitors <- NormalizeData(object = hypo_progenitors, normalization.method = "LogNormalize", scale.factor = 10000)
hypo_progenitors <- ScaleData(object = hypo_progenitors)
hypo_progenitors <- FindVariableGenes(object = hypo_progenitors)
length(x = hypo_progenitors@var.genes)

hypo_progenitors <- RunPCA(object = hypo_progenitors, pc.genes = hypo_progenitors@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)

pdf("Figures/PCElbowPlot_hypo_progenitors.pdf")
PCElbowPlot(object = hypo_progenitors, num.pc = 100) # Looks like anywhere from 20-30
dev.off()

hypo_progenitors <- FindClusters(object = hypo_progenitors, reduction.type = "pca", dims.use = 1:30, resolution = 1.2, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

hypo_progenitors <- RunTSNE(object = hypo_progenitors, dims.use = 1:30, check_duplicates = FALSE, do.fast = TRUE)

pdf("Figures/hypo_progenitors_cluster_plot.pdf",  height = 7, width = 7)
hypo_progenitors_plot <- TSNEPlot(object = hypo_progenitors, do.label = TRUE, no.legend = TRUE)
dev.off()

FeaturePlot(object = hypo_progenitors, features.plot = c("sox2"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 2)




png("Figures/hypo_progenitors_cluster_plot_gfap.png", height = 250, width = 250)
FeaturePlot(object = hypo_progenitors, features.plot = c("gfap"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, pch.use = 1)
dev.off()
png("Figures/hypo_progenitors_cluster_plot_sox2.png", height = 250, width = 250)
FeaturePlot(object = hypo_progenitors, features.plot = c("sox2"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, pch.use = 1)
dev.off()
png("Figures/hypo_progenitors_cluster_plot_fabp7a.png", height = 250, width = 250)
FeaturePlot(object = hypo_progenitors, features.plot = c("fabp7a"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, pch.use = 1)
dev.off()
png("Figures/hypo_progenitors_cluster_plot_cdh2.png", height = 250, width = 250)
FeaturePlot(object = hypo_progenitors, features.plot = c("cdh2"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, pch.use = 1)
dev.off()
png("Figures/hypo_progenitors_cluster_plot_tnc.png", height = 250, width = 250)
FeaturePlot(object = hypo_progenitors, features.plot = c("tnc"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, pch.use = 1)
dev.off()







# Find markers

# hypo.markers.10.12 <- FindMarkers(hypo_progenitors, ident.1 = 10, ident.2 = 12, min.pct = .25)

hypo_progenitors.markers <- FindAllMarkers(object = hypo_progenitors, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# Generate Feature plot of top 2 markers each

features.markers <- hypo_progenitors.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

pdf("Figures/cluster.markers.hypo_progenitors.pdf")
FeaturePlot(object = hypo_progenitors, features.plot = features.markers$gene, min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf("Figures/cluster.markers.vlnplot.hypo_progenitors.pdf")
VlnPlot(object = hypo_progenitors, features.plot = features.markers$gene, point.size.use = NA)
dev.off()

# Generate lists + feature plots of top 20 markers for each cluster

for (i in 0:12) {
  name <- paste("cluster_progenitors_", i, sep="")
  assign(name, hypo_progenitors.markers %>% filter(cluster == i))
}

# Generate feature plots for each cluster (top 20 genes)

for (i in 0:14) {
  cluster <- hypo_progenitors.markers %>% filter(cluster == i) %>% top_n(10, avg_logFC)
  pdf(paste("Figures/cluster_markers_progenitors/hypo.cluster.markers.hypo_progenitors.", i, ".pdf", sep = ""))
  FeaturePlot(object = hypo_progenitors, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.1)
  dev.off()
}


top10 <- hypo_progenitors.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top20 <- hypo_progenitors.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

pdf("Figures/cluster.markers.hypo_progenitors.heatmap.top10.pdf")
DoHeatmap(hypo_progenitors, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
pdf("Figures/cluster.markers.hypo_progenitors.heatmap.top20.pdf")
DoHeatmap(hypo_progenitors, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()