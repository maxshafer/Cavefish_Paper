library(Seurat)
library(Matrix)
library(dplyr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

# Load the datasets from 10x

data.vec <- c("raw_data/SF_M_2/AstMex102/", 
                 "raw_data/SF_M_3/AstMex102/",
                 "raw_data/SF_M_4/AstMex102/",
                 "raw_data/SF_M_5/AstMex102/",
                 "raw_data/SF_F_2/AstMex102/",
                 "raw_data/SF_F_3/AstMex102/",
                 "raw_data/SF_F_4/AstMex102/",
                 "raw_data/SF_F_5/AstMex102/",
                 "raw_data/MF_M_1/AstMex102/",
                 "raw_data/MF_F_1/AstMex102/",
                 "raw_data/TF_M_1/AstMex102/",
                 "raw_data/TF_F_1/AstMex102/",
                 "raw_data/PF_M_1/AstMex102/",
                 "raw_data/PF_F_1/AstMex102/",
                 "raw_data/PF_M_2/AstMex102/",
                 "raw_data/PF_F_2/AstMex102/")

names(data.vec) <- c("SFM2","SFM3","SFM4","SFM5","SFF2","SFF3","SFF4","SFF5","MFM1","MFF1","TFM1","TFF1","PFM1","PFF1","PFM2","PFF2")

input.data <- Read10X(data.vec)

# Create Seurat object, normalize, scale and find variable genes

hypo <- CreateSeuratObject(input.data, min.features = 200, project = "Astyanax_mexicanus") # 32191 genes across 66993 samples.

# Add meta data for sex, and population with correct identifiers
hypo@meta.data$species <- "astyanax_cave"
hypo@meta.data$species[grep("SF", hypo@meta.data$orig.ident)] <- "astyanax_surface"

hypo@meta.data$sex <- "male"
hypo@meta.data$sex[grep("FF", hypo@meta.data$orig.ident)] <- "female"

hypo@meta.data$morph <- "Choy_surface"
hypo@meta.data$morph[grep("PF", hypo@meta.data$orig.ident)] <- "Pachon_cave"
hypo@meta.data$morph[grep("TF", hypo@meta.data$orig.ident)] <- "Tinaja_cave"
hypo@meta.data$morph[grep("MF", hypo@meta.data$orig.ident)] <- "Molino_cave"


# Factor levels for meta data
hypo@meta.data$orig.ident <- factor(hypo@meta.data$orig.ident, levels = c("SFM2","SFM3","SFM4","SFM5","SFF2","SFF3","SFF4","SFF5","MFM1","MFF1","TFM1","TFF1","PFM1","PFF1","PFM2","PFF2"))
hypo@meta.data$sex <- factor(hypo@meta.data$sex, levels = c("female", "male"))


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
# 4 cells have the exact same PC values, so need to add check_duplicates = F to RunTSNE call - due to chance
# row.names(hypo@reductions$pca@cell.embeddings[duplicated(hypo@reductions$pca@cell.embeddings) | duplicated(hypo@reductions$pca@cell.embeddings, fromLast = TRUE),1:10])
# [1] "MFM1_CGAACATCAGCATGAG-1" "TFM1_TCAGCAAAGATAGTCA-1" "TFM1_TTTGCGCAGTTGTCGT-1" "PFF1_TACTTGTTCCGCTGTT-1"

hypo <- RunTSNE(object = hypo, reduction = "pca", dims = 1:50, tsne.method = "Rtsne", reduction.name = "tsne", reduction.key = "tsne_", seed.use = 1, check_duplicates = F)

hypo.ast <- hypo

save(hypo.ast, file = "AstMex_63k.Robj")
