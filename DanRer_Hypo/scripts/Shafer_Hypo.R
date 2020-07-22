library(Seurat)
library(Matrix)
library(dplyr)

setwd("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo")

# Load the datasets from 10x

data.vec <- c("raw_data/Hypo1_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo1_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo1_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo1_4/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_4/outs/filtered_gene_bc_matrices/dr86/",
              "raw_data/Hypo3_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo3_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo3_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo3_4/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_4/outs/filtered_gene_bc_matrices/dr86/")

names (data.vec) <- c("hypo1","hypo2","hypo3","hypo4","hypo5","hypo6","hypo7","hypo8","hypo9","hypo10","hypo11","hypo12","hypo13","hypo14","hypo15","hypo16")

input.data <- Read10X(data.vec)

# Create Seurat object, normalize, scale and find variable genes

hypo <- CreateSeuratObject(input.data, min.genes = 200, project = "Shafer_Hypo") # 32191 genes across 66993 samples.

# Add meta data for sex, and populate with correct identifiers

hypo@meta.data$sex <- "male"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo5"] <- "female1"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo6"] <- "female1"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo7"] <- "female1"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo8"] <- "female1"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo9"] <- "male2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo10"] <- "male2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo11"] <- "male2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo12"] <- "male2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo13"] <- "female2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo14"] <- "female2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo15"] <- "female2"
hypo@meta.data$sex[hypo@meta.data$orig.ident == "hypo16"] <- "female2"

# female_idents <- c("hypo5","hypo6","hypo7","hypo8","hypo13","hypo14","hypo15","hypo16")
# hypo_neuronal@meta.data$sex[hypo_neuronal@meta.data$orig.ident == female_idents] <- "female"

# QC figures

VlnPlot(hypo, features.plot = c("nGene", "nUMI"), group.by = "orig.ident", nCol = 2, point.size.use = 0.01)
VlnPlot(hypo, features.plot = c("nGene", "nUMI"), group.by = "sex", nCol = 2, point.size.use = 0.01)
GenePlot(hypo, gene1 = "nUMI", gene2 = "nGene")

# Filter out cells with nGene higher then 3k, nUMI higher then 15k?

hypo <- FilterCells(hypo, subset.names = c("nGene"), low.thresholds = c(200), high.thresholds = c(2500)) # 32191 genes across 65947 samples.

# Normalize, scale and find variable genes

hypo <- NormalizeData(object = hypo, normalization.method = "LogNormalize", scale.factor = 10000)
hypo <- FindVariableGenes(object = hypo, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = TRUE)
hypo <- ScaleData(object = hypo, genes.use = hypo@var.genes, do.par = TRUE, num.cores = 6)
length(x = hypo@var.genes)


# Run PCA for clustering

hypo <- RunPCA(object = hypo, pc.genes = hypo@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)

# # Determine the PCAs to use for clustering - using all 20
# 
# hypo <- JackStraw(hypo, num.pc = 150, num.replicate = 100, num.cores = 6)

pdf("Figures/PCElbowPlot.pdf", height = 14, width = 14)
PCElbowPlot(object = hypo, num.pc = 50)
dev.off()

pdf("Figures/PCHeatmap_1-5_45-50.pdf", height = 14, width = 14)
PCHeatmap(object = hypo, pc.use = c(1:5, 45:50), cellsuse = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()

pdf("Figures/PCHeatmap_1-5_70-75.pdf", height = 14, width = 14)
PCHeatmap(object = hypo, pc.use = c(1:5, 70:75), cellsuse = 500, do.balanced = TRUE, label.columns = FALSE)
dev.off()

# # pdf("Figures/JackStrawPlot_50.pdf")
# JackStrawPlot(hypo, PCs = 75:100)
# dev.off()


# Find clusters! Run for multiple resolutions

hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:50, resolution = 0.3, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
hypo <- FindClusters(object = hypo, resolution = 0.6)
hypo <- FindClusters(object = hypo, resolution = 1.2)
hypo <- FindClusters(object = hypo, resolution = 3.0)
hypo <- FindClusters(object = hypo, resolution = 5.0)
hypo <- FindClusters(object = hypo, resolution = 10.0)

# hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:75, resolution = 0.6, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
# hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:75, resolution = 1.2, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
# hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:75, resolution = 3.0, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
# hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:75, resolution = 5.0, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
# hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:75, resolution = 10.0, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)

# Calculate tSNE for plotting

hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, nthreads = 6, do.fast = TRUE)

set.seed(1188)
hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, tsne.method = "FIt-SNE", nthreads = 6, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "~/Documents/Programs/FIt-SNE-master/bin/fast_tsne", max_iter = 2000)

# cols <- c("lightgoldenrod1", "springgreen4")
cols <- c("firebrick1", "indianred1", "royalblue1", "skyblue1")
cols2 <- c("royalblue1", "skyblue2", "skyblue3", "skyblue4", "firebrick1", "firebrick2", "firebrick3", "firebrick4", "royalblue2", "royalblue3", "royalblue4", "indianred1", "indianred2", "indianred3", "indianred"," skyblue1")

# cluster <- DimPlot(hypo, reduction.use = "FItSNE", do.label = TRUE, no.legend = TRUE, group.by = "res.0.3", no.axes = TRUE, do.return = TRUE, vector.friendly = TRUE, pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
# sample <- DimPlot(hypo, reduction.use = "FItSNE", no.legend = FALSE, no.axes = TRUE, group.by = "orig.ident", do.return = TRUE, vector.friendly = TRUE, pt.size = 0.5, cols.use = cols2) + ggtitle("Sample") + theme(plot.title = element_text(hjust = 0.5))
# sex <- DimPlot(hypo, reduction.use = "FItSNE", no.legend = FALSE, no.axes = TRUE, group.by = "sex", do.return = TRUE, vector.friendly = TRUE, pt.size = 0.5, cols.use = cols) + ggtitle("Sex") + theme(plot.title = element_text(hjust = 0.5))

# plot_grid(p1, p2, p3)

# pdf("Figures/hypo_cluster_plot_res.0.3.pdf", height = 7, width = 7)
# DimPlot(hypo, reduction.use = "FItSNE", do.label = TRUE, no.legend = TRUE, group.by = "res.0.3", do.return = TRUE, vector.friendly = TRUE, pt.size = 1) + ggtitle("Cluster ID - res 0.3") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# pdf("Figures/hypo_cluster_plot_res.0.6.pdf", height = 7, width = 7)
# DimPlot(hypo, reduction.use = "FItSNE", do.label = TRUE, no.legend = TRUE, group.by = "res.0.6", do.return = TRUE, vector.friendly = TRUE, pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()
# pdf("Figures/hypo_cluster_plot_res.1.2.pdf", height = 7, width = 7)
# DimPlot(hypo, reduction.use = "FItSNE", do.label = TRUE, no.legend = TRUE, group.by = "res.1.2", do.return = TRUE, vector.friendly = TRUE, pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

# pdf("Figures/hypo_cluster_plot_orig_ident.pdf", height = 7, width = 8)
# DimPlot(hypo, reduction.use = "FItSNE", no.legend = FALSE, group.by = "orig.ident", do.return = TRUE, vector.friendly = TRUE, pt.size = 0.5) + ggtitle("Sample") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

# pdf("Figures/hypo_cluster_plot_sex.pdf", height = 7, width = 8)
# DimPlot(hypo, reduction.use = "FItSNE", no.legend = FALSE, group.by = "sex", do.return = TRUE, vector.friendly = TRUE, pt.size = 0.5, cols.use = cols) + ggtitle("Sex") + theme(plot.title = element_text(hjust = 0.5))
# dev.off()

## Make TSNE graphs

png("Figures/hypo_cluster_plot_res.0.3.pdf", height = 7, width = 7)
hypo_plot_res.0.3 <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.0.3", no.legend = TRUE, pt.size = .1)
dev.off()
png("Figures/hypo_cluster_plot_res.0.6_annotations.png", height = 7, width = 7, res = 250, units = "in")
hypo_plot_res.0.4 <- TSNEPlot(object = hypo, no.axes = TRUE, do.label = TRUE, group.by = "annotations_0.6", no.legend = TRUE, pt.size = .5)
dev.off()
png("Figures/hypo_cluster_plot_res.0.6.pdf", height = 7, width = 7)
hypo_plot_res.0.6 <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.0.6", no.legend = TRUE, pt.size = .1)
dev.off()
pdf("Figures/hypo_cluster_plot_res.1.2.pdf", height = 7, width = 7)
hypo_plot_res.1.2 <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.1.2", no.legend = TRUE, pt.size = .1)
dev.off()
pdf("Figures/hypo_cluster_plot_res.3.pdf", height = 7, width = 7)
hypo_plot_res.3 <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.3", no.legend = TRUE, pt.size = .1)
dev.off()
pdf("Figures/hypo_cluster_plot_res.5.pdf", height = 7, width = 7)
hypo_plot_res.5 <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.5", no.legend = TRUE, pt.size = .1)
dev.off()
pdf("Figures/hypo_cluster_plot_res.10.pdf", height = 7, width = 7)
hypo_plot_res.10 <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.10", no.legend = TRUE, pt.size = .1)
dev.off()

png("Figures/hypo_cluster_plot_subcluster.png", height = 8, width = 8, units = "in", res = 250)
hypo_plot_orig_ident <- TSNEPlot(object = hypo, group.by = "subcluster.NEW", pt.size = .05, do.label = T, no.legend = TRUE, no.axes = T)
dev.off()
png("Figures/hypo_cluster_plot_orig_ident.png", height = 8, width = 8, units = "in", res = 250)
hypo_plot_orig_ident <- TSNEPlot(object = hypo, group.by = "orig.ident", pt.size = .05, do.label = FALSE, no.legend = TRUE, no.axes = T, colors.use = cols2)
dev.off()
png("Figures/hypo_cluster_plot_sex.png", height = 8, width = 8, units = "in", res = 250)
hypo_plot_sex <- TSNEPlot(object = hypo, group.by = "sex", do.label = FALSE, pt.size = .5, no.legend = TRUE, no.axes = T, colors.use = cols)
dev.off()
png("Figures/hypo_cluster_fplot_vax1.png", height = 8, width = 8, units = "in", res = 250)
FeaturePlot(object = hypo, features.plot = c("vax1"), min.cutoff = "q9", cols.use = c("grey85", "blue", "blue", "blue"), no.axes = T, pt.size = .4)
dev.off()


FeaturePlot(object = hypo, no.axes = T, features.plot = c("hsp90aa1.2"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = 1)

DotPlot(object = hypo, genes.plot = c("lhx9", "NP5", "rfx4", "npvf", "hcrt"), plot.legend = TRUE, x.lab.rot = TRUE) + ylab(element_text(size = 2))


pdf("Figures/hypo_cluster_vlnplot_gaba_glut_res.0.3.pdf", height = 3, width = 12)
VlnPlot(hypo, features.plot = c("slc17a6a", "slc32a1"), nCol = 1, point.size.use = NA, size.title.use = 18)
dev.off()

pdf("Figures/hypo_cluster_vlnplot_neuronalres.0.3.pdf", height = 3, width = 12)
VlnPlot(hypo, features.plot = c("snap25a", "snap25b"), nCol = 1, point.size.use = NA)
dev.off()

pdf("Figures/hypo_cluster_vlnplot_neuronalres.0.3.pdf", height = 3, width = 12)
VlnPlot(hypo, features.plot = c("snap25a", "syt1a", "slc17a6a", "slc17a6b", "slc32a1", "th2", "olig2", "vim", "otpa", "prdx1", "her4.2", "adcyap1b"), nCol = 1, point.size.use = NA)
dev.off()

# # Computer overlap between var genes and CNE and duplicated genes from cichlid databases
# 
# pdf("Figures/hypo_cluster_vlnplot_cichlid_cne.pdf", height = 15, width = 12)
# VlnPlot(hypo, features.plot = c("cd151l", "lhx1a", "gsc", "epas1a", "uncx4.1", "prr15la", "epas1b", "foxp4", "nr5a1a", "scrt1a"), group.by = "res.3", nCol = 1, point.size.use = NA)
# dev.off()
# 
# pdf("Figures/hypo_cluster_vlnplot_cichlid_dup.pdf", height = 7.5, width = 12)
# VlnPlot(hypo, features.plot = c("grna", "ch25h", "cyp4t8", "si:dkey-238o13.4", "epdl1"), group.by = "res.5", nCol = 1, point.size.use = NA)
# dev.off()

# VlnPlot(hypo, features.plot = c("thy1", "ecscr"), group.by = "sex", nCol = 1, point.size.use = 1)

FeaturePlot(object = hypo, features.plot = c("hcrt"), cols.use = c("grey85", "blue", "blue", "blue"), min.cutoff = "q9", pt.size = .5, pch.use = 1)
FeaturePlot(object = hypo, features.plot = c("aspm", "fezf2", "rx3", "vim"), cols.use = c("grey85", "blue", "blue", "blue"), min.cutoff = "q9", pt.size = 1)
FeaturePlot(object = hypo, features.plot = c("hmx2", "hmx3a", "galn", "gsx1"), cols.use = c("grey85", "blue", "blue", "blue"), min.cutoff = "q9", pt.size = 1)
FeaturePlot(object = hypo, features.plot = c("crhb", "kiss2", "agrp", "npvf"), cols.use = c("grey85", "blue", "blue", "blue"), min.cutoff = "q9", pt.size = 1)

FeaturePlot(object = hypo, features.plot = c("snap25a", "snap25b"), cols.use = c("grey", "blue", "blue", "blue"), overlay = TRUE, no.legend = TRUE)
FeaturePlot(object = hypo, features.plot = c("sall3a", "foxo1a"), cols.use = c("grey85", "red", "blue", "green"), no.axes = T, overlay = TRUE, no.legend = TRUE, pt.size = .3)

# table_freq <- table(hypo@ident, hypo@meta.data$morph)
# freq_table_both <- table_freq[1:27,]
# freq_table_both <- prop.table(x = freq_table_both, margin = 2)
# freq_table_both <- freq_table_both[order(freq_table_both[,1]),]

freq_table <- prop.table(table(hypo@ident, hypo@meta.data$sex), margin = 2)

pdf("Figures/hypo_cluster_sex_freq_both.pdf", height = 30, width = 30)
barplot(height = freq_table, legend.text = T)
dev.off()

freq_table <- prop.table(x = table(hypo@meta.data$sex, hypo@ident), margin = 2)
freq_table <- freq_table[,order(freq_table[1,])]

pdf("Figures/hypo_cluster_sex_freq.pdf", height = 30, width = 30)
barplot(height = freq_table, legend.text = T)
dev.off()

# Find markers

# diff.markers <- FindMarkers(hypo, ident.1 = "43.0", ident.2 = "5.1", logfc.threshold = 0.25, min.pct = 0.25, only.pos = FALSE)
# dim(diff.markers)
# diff.markers[1:10,]

# Finding markers for cell clusters - small cluster fail to work against the whole dataset, so subset to 1k or 500 cells per cluster first

hypo <- SetAllIdent(hypo, id = "subcluster.NEW")
num.clusters <- (dim(table(hypo@ident)) - 1)

hypo.small <- SubsetData(hypo, max.cells.per.ident = 100)

# For single cluster
cluster.markers <- FindMarkers(hypo.small, ident.1 = 43, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.5, max.cells.per.ident = 500)

# For all clusters
# hypo.markers <- FindAllMarkers(object = hypo.small, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.5, max.cells.per.ident = 500)
hypo.markers.NEW <- FindAllMarkers(object = hypo.small, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.5)

write.csv(hypo.markers.NEW, file = "CSV/full_dataset/hypo.markers.res.0.4.csv")

# Generate Feature plot of top 2 markers each

features.markers <- hypo.markers.NEW %>% group_by(cluster) %>% top_n(1, avg_logFC)

pdf("Figures/hypo_cluster.markers.0.4.pdf", height = 60, width = 30)
FeaturePlot(object = hypo, features.plot = features.markers$gene, min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf("Figures/hypo_cluster.markers.0.4.vlnplot.pdf", height = 80, width = 20)
VlnPlot(object = hypo, features.plot = features.markers$gene, nCol = 1, size.x.use = 0, point.size.use = 0.01, size.title.use = 0)
dev.off()

pdf("Figures/hypo_cluster.markers.0.4.dotplot.pdf", height = 20, width = 20)
DotPlot(object = hypo, genes.plot = unique(features.markers$gene), plot.legend = TRUE, x.lab.rot = TRUE)
dev.off()

# Generate lists + feature plots of top 10 markers for each cluster

for (i in 0:num.clusters) {
  # name <- as.character(paste("cluster_", i, sep=""))
  hm <- filter(hypo.markers, cluster == i)
  # assign(name, hm)
  write.csv(hm, file = paste("CSV/full_dataset/cluster_markers/res.0.6.cluster_", i, ".csv", sep=""))
}

# Generate feature plots for each cluster (top 20 genes)

for (i in 0:43) {
  cluster <- hypo.markers %>% filter(cluster == i) %>% top_n(9, avg_logFC)
  pdf(paste("Figures/full_dataset/cluster_markers/hypo.cluster.markers.0.6.", i, ".pdf", sep = ""), height = 15, width = 15)
  FeaturePlot(object = hypo, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), nCol = 3, pt.size = 0.1)
  dev.off()
}


top10 <- hypo.markers.NEW %>% group_by(cluster) %>% top_n(10, avg_logFC)
top20 <- hypo.markers.NEW %>% group_by(cluster) %>% top_n(20, avg_logFC)

pdf("Figures/full_dataset/cluster.markers.0.4.heatmap.top10.pdf", height = 30, width = 100)
DoHeatmap(object = SubsetData(hypo, max.cells.per.ident = 50), genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
pdf("Figures/full_dataset/cluster.markers.0.4.heatmap.top20.pdf", height = 60, width = 100)
DoHeatmap(object = SubsetData(hypo, max.cells.per.ident = 50), genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()


# save seurat object and cluster markers

write.table(hypo.markers, file = "CSV/full_dataset/hypo.cluster.markers.0.6.txt", sep = "\t")



## Add annotation id labels to seurat object

test <- read.csv("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/CSV/full_dataset/Cluster_annotations_Shafer_Hypo.res.0.6.csv", header = TRUE, row.names = "cluster")

current.cluster.ids <- c(0:43)
new.cluster.ids <- as.vector(test$Subtype)
# new.cluster.ids <- as.vector(test[,1])

hypo <- SetAllIdent(hypo, id = "res.0.6")
hypo@ident <- mapvalues(x = hypo@ident, from = current.cluster.ids, to = new.cluster.ids)

hypo <- StashIdent(hypo, save.name = "annotations_0.6")
# hypo <- StashIdent(hypo, save.name = "neuronal")


save(hypo, file = "Shafer_Hypo_66k.Robj")

############################################################################
################### Subset seurat object for just neuronal or non-neuronal
############################################################################

neuronal_cells <- WhichCells(hypo, ident = c(0, 2:7, 10:12, 19:20, 22:24))

progenitor_cells <- WhichCells(hypo, ident = c(4))

nonneuronal_cells <- WhichCells(hypo, ident = c(8, 13:18))

input.data.neuronal <- hypo@raw.data[,neuronal_cells]
input.data.progenitors <- hypo@raw.data[,progenitor_cells]
input.data.nonneuronal <- hypo@raw.data[,nonneuronal_cells]

# Generate matrices and graphs for neuronal cells only
source("Neuronal_cells.R")

# Generate matrices and graphs for progenitor cells only
source("Progenitor_cells.R")

# Generate matrices and graphs for non-neuronal cells only
source("Nonneuronal_cells.R")


# Isolate only NP5 cells

hcrt.cells <- WhichCells(hypo.zeb, subset.name = "hcrt", accept.low = 1)
hypo.zeb <- SetIdent(hypo.zeb, cells.use = hcrt.cells, ident.use = 45)
hcrt.markers <- FindMarkers(hypo.zeb, ident.1 = "male", ident.2 = 45, max.cells.per.ident = 2000)
hcrt <- SubsetData(hypo.zeb, cells.use = hcrt.cells, subset.raw = T)

hcrt.raw <- as.matrix(hcrt.raw)
write.csv(NP5.raw, file = "hcrt_raw_data.csv")
write.csv(NP5.markers, file = "hcrt_diff_exp.csv")


list <- read.csv("test.csv", header = TRUE, stringsAsFactors = FALSE)

hypo@meta.data$subcluster.NEW <- list$NEW[match(hypo@meta.data$subcluster.0.4, list$OLD)]



