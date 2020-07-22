library(Seurat)
library(Matrix)
library(dplyr)

setwd("~/Documents/Schier_Lab/R_Projects/AstMex_Hypo")

# Load the datasets from 10x

load("AstMex_cf.Robj")
load("AstMex_sf.Robj")

# Create Seurat object, normalize, scale and find variable genes

hypo <- MergeSeurat(object1 = hypo.sf, object2 = hypo.cf, do.normalize = FALSE)
hypo <- NormalizeData(hypo, normalization.method = "TFIDF")
# hypo <- FindVariableGenes(object = hypo, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = TRUE)
hypo@var.genes <- FindTopTFIDF(hypo, gene.number = 600)
hypo <- ScaleData(object = hypo, genes.use = hypo@var.genes, do.par = TRUE, num.cores = 6)
length(x = hypo@var.genes)

# QC figures

VlnPlot(hypo, features.plot = c("nGene", "nUMI"), group.by = "orig.ident", nCol = 2, point.size.use = 0.01)
VlnPlot(hypo, features.plot = c("nGene", "nUMI"), group.by = "sex", nCol = 2, point.size.use = 0.01)
GenePlot(hypo, gene1 = "nUMI", gene2 = "nGene")


# Run PCA for clustering

hypo <- RunPCA(object = hypo, pc.genes = hypo@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)

# # Determine the PCAs to use for clustering - using all 20
# 
# hypo <- JackStraw(hypo, num.pc = 150, num.replicate = 100, num.cores = 6)

# pdf("Figures/PCElbowPlot.pdf", height = 14, width = 14)
PCElbowPlot(object = hypo, num.pc = 150)
# dev.off()

# pdf("Figures/PCHeatmap_1-5_45-50.pdf", height = 14, width = 14)
# PCHeatmap(object = hypo, pc.use = c(1:5, 45:50), cellsuse = 500, do.balanced = TRUE, label.columns = FALSE)
# dev.off()

# pdf("Figures/PCHeatmap_1-5_70-75.pdf", height = 14, width = 14)
# PCHeatmap(object = hypo, pc.use = c(1:5, 70:75), cellsuse = 500, do.balanced = TRUE, label.columns = FALSE)
# dev.off()

# # pdf("Figures/JackStrawPlot_50.pdf")
# JackStrawPlot(hypo, PCs = 75:100)
# dev.off()


# Find clusters! Run for multiple resolutions

hypo <- FindClusters(object = hypo, reduction.type = "pca", dims.use = 1:50, resolution = 0.3, save.SNN = TRUE, n.start = 10, nn.eps = 0.5, print.output = FALSE, force.recalc = TRUE)
hypo <- FindClusters(object = hypo, resolution = 0.6)
hypo <- FindClusters(object = hypo, resolution = 1.2)
hypo <- FindClusters(object = hypo, resolution = 3.0)
# hypo <- FindClusters(object = hypo, resolution = 5.0)
# hypo <- FindClusters(object = hypo, resolution = 10.0)


# Calculate tSNE for plotting

# hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, nthreads = 4, do.fast = TRUE, check_duplicates = FALSE)

set.seed(1188)
hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, tsne.method = "FIt-SNE", nthreads = 6, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "~/Downloads/FIt-SNE-ec25f1b36598a2d21869d10a258ac366a12f0b05/bin/fast_tsne", max_iter = 2000)

hypo.ast <- hypo

save(hypo.ast, file = "AstMex_66k.Robj")