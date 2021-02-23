library(Seurat)
library(patchwork)
library(ggplot2)

# Load the initial Seurat objects

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

# Load the datasets from 10x

data.vec.zeb <- c("raw_data/Hypo1_1/outs/filtered_gene_bc_matrices/dr86/", 
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

names (data.vec.zeb) <- c("hypo1","hypo2","hypo3","hypo4","hypo5","hypo6","hypo7","hypo8","hypo9","hypo10","hypo11","hypo12","hypo13","hypo14","hypo15","hypo16")
input.data <- Read10X(data.vec.zeb)

# Create Seurat object, normalize, scale and find variable genes

hypo.zeb.og <- CreateSeuratObject(input.data, min.features = 200, project = "Danio_rerio")


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

# Load the datasets from 10x

data.vec.ast <- c("raw_data/SF_M_2/AstMex102/", 
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

names(data.vec.ast) <- c("SFM2","SFM3","SFM4","SFM5","SFF2","SFF3","SFF4","SFF5","MFM1","MFF1","TFM1","TFF1","PFM1","PFF1","PFM2","PFF2")

input.data <- Read10X(data.vec.ast)

# Create Seurat object, normalize, scale and find variable genes

hypo.ast.og <- CreateSeuratObject(input.data, min.features = 200, project = "Astyanax_mexicanus")

# Save files for loading
saveRDS(hypo.zeb.og, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_og.rds")
saveRDS(hypo.ast.og, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_og.rds")

hypo.zeb.og <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_og.rds")
hypo.ast.og <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_og.rds")
## Run QC + make figures

hypo.zeb.og[["percent.mt"]] <- PercentageFeatureSet(hypo.zeb.og, pattern = "^mt-")
hypo.ast.og[["percent.mt"]] <- PercentageFeatureSet(hypo.ast.og, pattern = "^mt-")

## Make final plots

new_theme <- function() {
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_blank(), 
    legend.position = "none"
  )
}

vln.plots <- VlnPlot(hypo.zeb.og, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of genes") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

vln.plots[[1]] <- VlnPlot(hypo.zeb.og, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of genes") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln.plots[[2]] <- VlnPlot(hypo.zeb.og, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of counts") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln.plots[[3]] <- VlnPlot(hypo.zeb.og, features = c("percent.mt"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("% mito.")

vln.plots[[4]] <- VlnPlot(hypo.ast.og, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of genes") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln.plots[[5]] <- VlnPlot(hypo.ast.og, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme()+ ylab("# of counts")

zeb.fs1 <- FeatureScatter(hypo.zeb.og, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.05) + new_theme() + theme(legend.position = c(0.85,0.35))
zeb.fs2 <- FeatureScatter(hypo.zeb.og, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.05) + new_theme() + theme(legend.position = c(0.75,0.55))

ast.fs1 <- FeatureScatter(hypo.ast.og, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.05) + new_theme() + theme(legend.position = c(0.85,0.35))

# Plot
layout <- "
AAAAAA
BBBBBB
CCCCCC
DDDDDD
EEEEEE
FFGGHH
FFGGHH
"
dev.new()
vln.plots + plot_layout(design = layout, guides = "collect", widths = unit(c(15), "mm"), height = unit(c(15), "mm"))

## Need to copy to affinity and rasterize one at a time (60k dots each)
layout2 <- "
FFGGHH
FFGGHH
"
dev.new()
zeb.fs1 + zeb.fs2 + ast.fs1 + plot_layout(design = layout2, guides = "collect", widths = unit(c(15), "mm"), height = unit(c(15), "mm"))
