library(Seurat)
library(patchwork)
library(networkD3)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")
Idents(hypo.integrated) <- "integrated_Cluster"

######### Non-homologous Genes #########

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_orthologs_corrected.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_orthologs.txt", head = TRUE)
names(mart) <- c("danio", "astyanax")

## Extract genes don't have a homolog in the other species
ast.genes <- mart[[2]][mart[[2]]$Zebrafish.gene.stable.ID == "",]
ast.genes$Gene.name[ast.genes$Gene.name == ""] <- ast.genes$Gene.stable.ID[ast.genes$Gene.name == ""]
ast.genes <- ast.genes$Gene.name

zeb.genes0 <- mart[[1]][mart[[1]]$Cave.fish.gene.stable.ID == "",]
zeb.genes1 <- zeb.genes0
zeb.genes1$Gene.name[!(is.na(zeb.genes1$GeneID))] <- zeb.genes1$GeneID[!(is.na(zeb.genes1$GeneID))]
zeb.genes <- zeb.genes1$Gene.name

# Some index files
clusters <- c("Neuronal_00", "Neuronal_01", "Neuronal_02", "Neuronal_03", "Neuronal_04", "Neuronal_05", "Neuronal_06", "Neuronal_07", "Neuronal_08", "Neuronal_09", "Neuronal_10", "Neuronal_11", "Neuronal_12", "Neuronal_13")
int.idents <- c(50, 25, 25, 25, 25, 20, 20, 15, 15, 15, 10, 10, 5, 5)

idents <- c("Neuronal_03","Neuronal_05", "Neuronal_07", "Neuronal_10", "Neuronal_12", "Neuronal_13")

subsets <- list()
# subsets.markers <- list()
for (i in 1:length(idents)) {
  subset <- subset(hypo.integrated, idents = idents[[i]])
  Idents(subset) <- "species.2"
  reference.list <- c(subset(subset, idents = "zebrafish"), subset(subset, idents = "astyanax"))
  for (j in 1:length(reference.list)) {
    reference.list[[j]] <- NormalizeData(object = reference.list[[j]], verbose = FALSE)
    reference.list[[j]] <- FindVariableFeatures(object = reference.list[[j]], selection.method = "vst", nfeatures = 1500, verbose = FALSE)
  }
  
  VariableFeatures(reference.list[[1]]) <- VariableFeatures(reference.list[[1]])[!(VariableFeatures(reference.list[[1]]) %in% zeb.genes)]
  VariableFeatures(reference.list[[2]]) <- VariableFeatures(reference.list[[2]])[!(VariableFeatures(reference.list[[2]]) %in% ast.genes)]
  
  # Integrated subsets
  subset.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:int.idents[match(idents[[i]], clusters)])
  subset.integrated <- IntegrateData(anchorset = subset.anchors, dims = 1:int.idents[match(idents[[i]], clusters)], k.weight = min(100, min(ncol(reference.list[[1]]@assays$RNA), ncol(reference.list[[2]]@assays$RNA))))
  
  
  # Scale and run dimensionality reduction on integrated data
  DefaultAssay(object = subset.integrated) <- "integrated"
  subset.integrated <- ScaleData(object = subset.integrated, verbose = FALSE)
  subset.integrated <- RunPCA(object = subset.integrated, npcs = int.idents[match(idents[[i]], clusters)], verbose = FALSE)
  subset.integrated <- RunTSNE(object = subset.integrated, reduction = "pca", dims = 1:int.idents[match(idents[[i]], clusters)], check_duplicates = F)
  
  subset.integrated <- FindNeighbors(subset.integrated, reduction = "pca", dims = 1:length(subset.integrated@reductions$pca))
  subset.integrated <- FindClusters(subset.integrated, resolution = 0.2, print.output = 0)
  subset.integrated <- FindClusters(subset.integrated, resolution = 0.25, print.output = 0)
  subset.integrated <- FindClusters(subset.integrated, resolution = 0.4, print.output = 0)
  
  # DefaultAssay(subset.integrated) <- "RNA"
  # subset.integrated.markers <- FindAllMarkers(subset.integrated)
  
  subsets[[i]] <- subset.integrated
  # subsets.markers[[i]] <- subset.integrated.markers
}


plots <- list()
for (i in 1:length(subsets)) {
  p1 <- DimPlot(subsets[[i]], reduction = "tsne", group.by = "integrated_Subcluster", label = T, label.size = 2, pt.size = 0.1, repel = T) + NoLegend() + NoAxes()
  p3 <- DimPlot(subsets[[i]], reduction = "tsne", group.by = "integrated_snn_res.0.25", label = T, label.size = 2, pt.size = 0.1, repel = T) + NoLegend() + NoAxes()
  plots[[i]] <- p1
  plots[[i+8]] <- p3
}
plots <- lapply(plots, function(x) x & theme(axis.text = element_blank(), axis.title = element_blank()))
# plots[[1]] <- plots[[1]] & theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10))

patchwork <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plots[[5]] + plots[[6]] + plots[[7]] + plots[[8]] + plots[[9]] + plots[[10]] + plots[[11]] + plots[[12]] + plots[[13]] + plots[[14]] + plots[[15]] + plots[[16]]
patchwork <- patchwork + plot_layout(byrow = F, ncol = 2, height = unit(c(22), c("mm")), width = unit(c(22), c("mm")), guides = "collect")

dev.new()
patchwork

# Make Sankey for comparing old subclusters, to new subclusters
data <- lapply(subsets, function(x) reshape2::melt(table(x@meta.data$integrated_Subcluster, x@meta.data$integrated_snn_res.0.25)))

nodes <- lapply(data, function(x) unique(c(as.character(unique(x$Var2)), as.character(unique(x$Var1)))))
nodes <- lapply(nodes, function(x) as.data.frame(x))

for (i in 1:length(data)) {
  data[[i]]$Var1 <- as.numeric(row.names(nodes[[i]])[match(data[[i]]$Var1, nodes[[i]]$x)])-1
  data[[i]]$Var2 <- as.numeric(row.names(nodes[[i]])[match(data[[i]]$Var2, nodes[[i]]$x)])-1
}

data2 <- lapply(data, function(x) x[x$value > 5,])

lapply(seq_along(data2), function(x) sankeyNetwork(Links = data2[[x]], Nodes = nodes[[x]], Source = "Var1", Target = "Var2", Value = "value", NodeID = "x", units = "cells", fontSize = 8, nodeWidth = 40, nodePadding = 5, height = 175, width = 350, iterations = 10000, sinksRight = F))


sankeyNetwork(Links = data[[2]], Nodes = nodes[[2]], Source = "Var1", Target = "Var2", Value = "value", NodeID = "x", units = "cells", fontSize = 8, nodeWidth = 40, nodePadding = 10, height = 250, width = 350, iterations = 10000, sinksRight = F)

#sankeyNetwork(Links = data2[[5]], Nodes = nodes[[5]], Source = "Var1", Target = "Var2", Value = "value", NodeID = "x", units = "cells", fontSize = 14)


