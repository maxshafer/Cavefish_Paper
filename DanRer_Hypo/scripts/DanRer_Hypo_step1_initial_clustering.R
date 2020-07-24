library(Seurat)
library(Matrix)
library(dplyr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

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

hypo <- CreateSeuratObject(input.data, min.features = 200, project = "Danio_rerio") # 32191 genes across 66993 samples.

# Add meta data for sex, and populate with correct identifiers

hypo@meta.data$sex <- "male1"
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

hypo@meta.data$species <- "zebrafish"


# Factor levels for meta data
hypo@meta.data$orig.ident <- factor(hypo@meta.data$orig.ident, levels = c("hypo5","hypo6","hypo7","hypo8","hypo13","hypo14","hypo15","hypo16","hypo1","hypo2","hypo3","hypo4","hypo9","hypo10","hypo11","hypo12"))
hypo@meta.data$sex <- factor(hypo@meta.data$sex, levels = c("female1", "female2", "male1", "male2"))


# Create Seurat object, normalize, scale and find variable genes

hypo <- NormalizeData(hypo, normalization.method = "LogNormalize", scale.factor = 10000)
hypo <- FindVariableFeatures(hypo, selection.method = "mvp")
length(VariableFeatures(hypo))
# [1] 488

hypo <- ScaleData(object = hypo, features = VariableFeatures(hypo))


# QC figures and subsetting

VlnPlot(hypo, features = c("nFeature_RNA", "nCount_RNA"), group.by = "orig.ident", ncol = 2)
FeatureScatter(hypo, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

hypo <- subset(hypo, subset = nFeature_RNA < 2500)

# Run PCA for clustering

hypo <- RunPCA(object = hypo, features = VariableFeatures(hypo), npcs = 100, set.seed = 0) # Reproduces! Needed to keep all the cells, doh!

# # Determine the PCAs to use for clustering - using all 20
# 
# hypo <- JackStraw(hypo, num.pc = 150, num.replicate = 100, num.cores = 6)

# pdf("Figures/PCElbowPlot.pdf", height = 14, width = 14)
# PCElbowPlot(object = hypo, num.pc = 150)
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


# Find clusters! Run using 50 PCs at resolution 0.6

hypo <- FindNeighbors(hypo, dims = 1:50, k.param = 30, nn.eps = 0.5)
hypo <- FindClusters(hypo, resolution = 0.6, random.seed = 0)

# Calculate tSNE for plotting

hypo <- RunTSNE(object = hypo, reduction = "pca", dims = 1:50, tsne.method = "Rtsne", reduction.name = "tsne", reduction.key = "tsne_", seed.use = 1, check_duplicates = F)


saveRDS(hypo, file = "DanRer_66k.rds")
