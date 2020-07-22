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


# Calculate tSNE for plotting

hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, nthreads = 6, do.fast = TRUE)

set.seed(1188)
hypo <- RunTSNE(object = hypo, reduction.use = "pca", dims.use = 1:50, tsne.method = "FIt-SNE", nthreads = 6, reduction.name = "FItSNE", reduction.key = "FItSNE_", fast_tsne_path = "~/Documents/Programs/FIt-SNE-master/bin/fast_tsne", max_iter = 2000)

save(hypo, file = "Shafer_Hypo_67k.Robj")




