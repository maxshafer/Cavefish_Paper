###################################### Neuronal cells ######################################

# hypo <- CreateSeuratObject(input.data.neuronal, min.genes = 200, project = "Shafer_Hypo")
hypo <- NormalizeData(object = hypo, normalization.method = "TFIDF")
# hypo <- FindVariableGenes(object = hypo, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = TRUE)
hypo@var.genes <- FindTopTFIDF(hypo, gene.number = 500)
length(x = hypo@var.genes)
hypo <- ScaleData(object = hypo, genes.use = hypo@var.genes, do.par = TRUE, num.cores = 6)

# Run PCA and find clusters

hypo <- RunPCA(object = hypo, pc.genes = hypo@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)

# Determine the PCAs to use for clustering - using all 20

hypo <- JackStraw(hypo, num.pc = 75, num.replicate = 75, do.par = T, num.cores = 6)
JackStrawPlot(hypo, PCs = 1:50)

pdf("Figures/PCElbowPlot_hypo.pdf")
PCElbowPlot(object = hypo, num.pc = 100) # Looks like anywhere from 20-30
dev.off()

hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:50, resolution = 0.3, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
hypo <- FindClusters(object = hypo, resolution = 0.6)
hypo <- FindClusters(object = hypo, resolution = 1.2)
hypo <- FindClusters(object = hypo, resolution = 3.0)
hypo <- FindClusters(object = hypo, resolution = 5.0)
hypo <- FindClusters(object = hypo, resolution = 10.0)

# Calculate tSNE for plotting

# hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, nthreads = 4, do.fast = TRUE, check_duplicates = FALSE)

set.seed(1188)
hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, tsne.method = "FIt-SNE", nthreads = 6, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "~/Downloads/FIt-SNE-ec25f1b36598a2d21869d10a258ac366a12f0b05/bin/fast_tsne", max_iter = 2000)

hypo <- RunUMAP(object = hypo, reduction.use = "pca", dims.use = 1:50, min_dist = .05)


# Make cluster plots

cluster <- DimPlot(hypo, group.by = "res.0.3", reduction.use = "FItSNE", do.label = TRUE, no.legend = TRUE, do.return = FALSE, vector.friendly = FALSE, pt.size = .5) + ggtitle("Cluster ID - res 0.3") + theme(plot.title = element_text(hjust = 0.5))


FeaturePlot(object = hypo, features.plot = c("slc17a6b", "slc17a6a", "slc32a1", "snap25b"), cols.use = c("lightgrey", "blue"), pt.size = 1)
FeaturePlot(object = hypo, reduction.use = "FItSNE", features.plot = c("gad1b"), cols.use = c("lightgrey", "blue"), pt.size = 1)











# VlnPlot(object = hypo, features.plot = c("NP5"), group.by = "res.3", point.size.use = 0.1)

gaba_clusters <- c(0, 2, 5:8, 13, 15:17, 22, 24, 26, 28, 30, 31, 33)
glut_clusters <- c(3, 6, 9, 10, 12, 20, 21, 23, 32)
other_clusters <- c(1, 4, 11, 14, 18, 19, 25, 27, 29)

hypo_gaba_neurons <- SubsetData(hypo, ident.use = gaba_clusters, subset.raw = T)
hypo_glut_neurons <- SubsetData(hypo, ident.use = glut_clusters, subset.raw = T)
hypo_other_neurons <- SubsetData(hypo, ident.use = other_clusters, subset.raw = T)

#### Find markers

hypo.markers <- FindAllMarkers(object = hypo, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# Generate Feature plot of top 2 markers each

features.markers <- hypo.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

pdf("Figures/cluster.markers_hypo.pdf")
FeaturePlot(object = hypo, features.plot = features.markers$gene, min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf("Figures/cluster.markers.vlnplot_hypo.pdf")
VlnPlot(object = hypo, features.plot = features.markers$gene, point.size.use = NA)
dev.off()

# Generate lists + feature plots of top 20 markers for each cluster

for (i in 0:25) {
  name <- paste("cluster_neuronal_", i, sep="")
  assign(name, hypo.markers %>% filter(cluster == i))
}

# Generate feature plots for each cluster (top 10 genes)

for (i in 0:25) {
  cluster <- hypo.markers %>% filter(cluster == i) %>% top_n(10, avg_logFC)
  pdf(paste("Figures/hypo.cluster.markers.hypo.", i, ".pdf", sep = ""))
  FeaturePlot(object = hypo, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.1)
  dev.off()
}


top10 <- hypo.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
top20 <- hypo.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)

pdf("Figures/cluster.markers.hypo.heatmap.top10.pdf")
DoHeatmap(hypo, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
pdf("Figures/cluster.markers.hypo.heatmap.top20.pdf")
DoHeatmap(hypo, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()