library(Seurat)
library(Matrix)
library(dplyr)

setwd("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo")


# Load the datasets from 10x
data.vec <- c("raw_data/hypo1_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/hypo1_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/hypo1_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/hypo1_4/outs/filtered_gene_bc_matrices/dr86/")

names (data.vec) <- c("hypo1", "hypo2", "hypo3", "hypo")

input.data <- Read10X(data.vec)

# Load Bushra's Hypothalamus data from Nature Biotech Paper

hypo.data <- as.matrix(read.table("raw_data/Max_ds.txt", row.names = 1, header = TRUE))

load("~/Documents/Schier_Lab/R_Projects/Raj_2017/Raj_seurat.Robj") # Called hypo :(

raj <- hypo
raj@meta.data$dataset <- "Raj"

# Create Seurat objects - min.genes = 200 gives

shafer <- CreateSeuratObject(input.data, min.genes = 200, project = "Shafer_Hypo")
shafer@meta.data$dataset <- "Shafer"

# raj <- CreateSeuratObject(raw.data = hypo.data, min.genes = 200, project = "Raj_Hypo)")
# raj@meta.data$dataset <- "Raj"


hypo <- MergeSeurat(shafer, raj, add.cell.id1 = "Shafer", add.cell.id2 = "Raj", project = "Hypo")


# hypo <- NormalizeData(object = hypo, normalization.method = "LogNormalize", scale.factor = 10000)
hypo <- ScaleData(object = hypo)
hypo <- FindVariableGenes(object = hypo)
length(x = hypo@var.genes)

write.csv(hypo@var.genes, file = "CSV/hypo.var.genes.csv")

hypo <- RunPCA(object = hypo, pc.genes = hypo@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)

# Determine the PCAs to use for clustering - using all 20

pdf("Figures/full_dataset/PCElbowPlot.pdf")
PCElbowPlot(object = hypo, num.pc = 100) # Looks like anywhere from 20-30
dev.off()


# Find clusters!

hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:50, resolution = 1.2, print.output = 0, save.SNN = TRUE)

hypo <- RunTSNE(object = hypo, dims.use = 1:50, check_duplicates = FALSE, do.fast = TRUE)

pdf("Figures/full_dataset/hypo_cluster_plot.pdf")
hypo_plot <- TSNEPlot(object = hypo, do.label = TRUE)
dev.off()

pdf("Figures/full_dataset/hypo_cluster_plot_dataset.pdf")
hypo_plot_dataset <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "dataset")
dev.off()

pdf("Figures/full_dataset/hypo_cluster_plot_orig.ident.pdf")
hypo_plot_orig.ident <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "orig.ident")
dev.off()

pdf("Figures/full_dataset/hypo_cluster_plot_raj.ident.pdf")
hypo_plot_raj.ident <- TSNEPlot(object = hypo, do.label = TRUE, group.by = "res.0.6")
dev.off()

VlnPlot(object = hypo, features.plot = c("snap25a", "snap25b", "syt1a"), point.size.use = NA)
VlnPlot(object = hypo, features.plot = c("snap25a", "snap25b", "syt1a"), point.size.use = NA)
FeaturePlot(object = hypo, features.plot = c("snap25a", "snap25b", "syt1a"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = .5)
FeaturePlot(object = hypo, features.plot = c("otpb", "prdx1", "prdx5"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = .5)

FeaturePlot(object = hypo, features.plot = c("prdx1", "npvf"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE)
FeaturePlot(object = hypo, features.plot = c("prdx1", "crhb"), cols.use = c("grey", "red", "blue", "green"), overlay = TRUE, no.legend = FALSE)








# # Find markers
# 
# hypo.markers.36 <- FindMarkers(hypo, ident.1 = 36, min.pct = .25)
# 
# hypo.markers <- FindAllMarkers(object = hypo, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
# 
# # Generate Feature plot of top 2 markers each
# 
# features.markers <- hypo.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)
# 
# pdf("Figures/full_dataset/cluster.markers.pdf")
# FeaturePlot(object = hypo, features.plot = features.markers$gene, min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
# dev.off()
# 
# pdf("Figures/full_dataset/cluster.markers.vlnplot.pdf")
# VlnPlot(object = hypo, features.plot = features.markers$gene, point.size.use = NA)
# dev.off()
# 
# # Generate lists + feature plots of top 10 markers for each cluster
# 
# for (i in 0:36) {
#   name <- paste("cluster_", i, sep="")
#   assign(name, hypo.markers %>% filter(cluster == i) %>% top_n(20, avg_logFC))
#   write.csv(cluster_markers, file = paste("CSV/full_dataset/cluster_markers/cluster_", i, ".csv", sep=""))
# }
# 
# # Generate feature plots for each cluster (top 20 genes)
# 
# for (i in 0:36) {
#   cluster <- hypo.markers %>% filter(cluster == i) %>% top_n(10, avg_logFC)
#   pdf(paste("Figures/full_dataset/cluster_markers/hypo.cluster.markers.", i, ".pdf", sep = ""))
#   FeaturePlot(object = hypo, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.1)
#   dev.off()
# }
# 
# 
# top10 <- hypo.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# top20 <- hypo.markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
# 
# pdf("Figures/full_dataset/cluster.markers.heatmap.top10.pdf")
# DoHeatmap(hypo, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
# dev.off()
# pdf("Figures/full_dataset/cluster.markers.heatmap.top20.pdf")
# DoHeatmap(hypo, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
# dev.off()



# # save seurat object and cluster markers

# write.table(hypo.markers, file = "cluster.markers.txt", sep = "\t")

# save(hypo, file = "Shafer_Hypo.Robj")

# ############################################################################
# ################### Subset seurat object for just neuronal or non-neuronal
# ############################################################################

# neuronal_cells <- WhichCells(hypo, ident = c(0, 1, 3, 4, 5, 6, 7, 8, 10, 12, 14, 15, 17, 21, 22, 23, 25, 28, 32, 34, 35, 36))

# progenitor_cells <- WhichCells(hypo, ident = c(2, 9, 16, 27, 31))

# nonneuronal_cells <- WhichCells(hypo, ident = c(11, 13, 18, 19, 20, 24, 26, 29, 30, 33))
  
# input.data.neuronal <- hypo@raw.data[,neuronal_cells]
# input.data.progenitors <- hypo@raw.data[,progenitor_cells]
# input.data.nonneuronal <- hypo@raw.data[,nonneuronal_cells]

# # Generate matrices and graphs for neuronal cells only
# source("Neuronal_cells.R")

# # Generate matrices and graphs for progenitor cells only
# source("Progenitor_cells.R")

# # Generate matrices and graphs for non-neuronal cells only
# source("Nonneuronal_cells.R")




