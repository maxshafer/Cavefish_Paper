library(randomForest)
library(Seurat)
library(Matrix)
library(dplyr)
library(reshape2)

plotConfusionMatrix = function(X, x.order=default.x.order, y.order=default.y.order, row.scale=TRUE, col.scale=FALSE, cols.use=gray.colors(10), max.size=5, ylab.use="Known", xlab.use="Predicted"){
  
  if (!col.scale & row.scale){ X = t(scale(t(X), center=FALSE, scale=rowSums(X)));  X=X*100 }
  if (col.scale & !row.scale){ X = scale(X, center=FALSE, scale=colSums(X)); X = X*100 }
  if(col.scale & row.scale){
    print("Only one of row.scale or col.scale should be true. performing row scaling by default")
    X = t(scale(t(X), center=FALSE, scale=rowSums(X)))
    X=X*100
  }
  default.x.order = colnames(X)
  default.y.order = rownames(X) # Default order is reverse (starts from axis 0)
  #X = X[rev(1:dim(X)[1]),]
  X = melt(X)
  colnames(X) = c("Known", "Predicted", "Percentage")
  X$Known = factor(X$Known, levels=rev(unique(X$Known)));
  X$Predicted = as.factor(X$Predicted)
  p = ggplot(X, aes(y = Known,  x = Predicted)) + geom_point(aes(colour = Percentage,  size =Percentage)) + 
    scale_color_gradient(low ="#edf8b1",   high = "#2c7fb8", limits=c(0, 100 ))+scale_size(range = c(1, max.size))+   theme_bw() #+nogrid
  p = p + xlab(xlab.use) + ylab(ylab.use) + theme(axis.text.x=element_text(size=12, face="italic", hjust=1)) + 
    theme(axis.text.y=element_text(size=12, face="italic"))
  p = p + scale_x_discrete(limit = x.order) + scale_y_discrete(limit = y.order)
  print(p)
}

setwd("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/30-04-2018_CCA")


# # Load the datasets from 10x
# data.vec <- c("raw_data/hypo1_1/outs/filtered_gene_bc_matrices/dr86/", 
#               "raw_data/hypo1_2/outs/filtered_gene_bc_matrices/dr86/", 
#               "raw_data/hypo1_3/outs/filtered_gene_bc_matrices/dr86/", 
#               "raw_data/hypo1_4/outs/filtered_gene_bc_matrices/dr86/")
# 
# names (data.vec) <- c("hypo1", "hypo2", "hypo3", "hypo")
# 
# input.data <- Read10X(data.vec)

# # Load Bushra's Hypothalamus data from Nature Biotech Paper
# 
# raj.data <- as.matrix(read.table("raw_data/Max_ds.txt", row.names = 1, header = TRUE))

# Load old objects and/or create new seurat objects - min.genes = 200 gives

hypo.original.pca <- hypo # backup my PCA clustered object
shafer <- hypo 

# shafer <- CreateSeuratObject(input.data, min.genes = 200, project = "Shafer_Hypo")
# shafer@meta.data$dataset <- "Shafer"
# shafer <- NormalizeData(object = shafer, normalization.method = "LogNormalize", scale.factor = 10000)
# shafer <- ScaleData(object = shafer)
# shafer <- FindVariableGenes(object = shafer)
# length(x = shafer@var.genes)

# raj <- CreateSeuratObject(raw.data = raj.data, min.genes = 200, project = "Raj_Hypo)")
# raj@meta.data$dataset <- "Raj"
# raj <- NormalizeData(object = raj, normalization.method = "LogNormalize", scale.factor = 10000)
# raj <- ScaleData(object = raj)
# raj <- FindVariableGenes(object = raj)
# length(x = raj@var.genes)

load("~/Documents/Schier_Lab/R_Projects/Raj_2017/Raj_seurat.Robj") # named hypo :(

raj <- hypo # rename Bushra's data

raj@meta.data$dataset <- "Raj"
hypo.original.pca@meta.data$dataset <- "Shafer"

var.union <- union(hypo.original.pca@var.genes, raj@var.genes)

hypo <- RunCCA(object = hypo.original.pca, object2 = raj, genes.use = var.union, num.cc = 50)

rm(raj)

# hypo <- MergeSeurat(shafer, raj, add.cell.id1 = "Shafer", add.cell.id2 = "Raj", project = "Hypo")

# visualize results of CCA plot CC1 versus CC2 and look at a violin plot
p <- DimPlot(object = hypo, reduction.use = "cca", group.by = "dataset", pt.size = 0.5, 
              do.return = TRUE)
p1 <- VlnPlot(object = hypo, features.plot = "CC1", group.by = "dataset", do.return = TRUE, point.size.use = 0)
p2 <- VlnPlot(object = hypo, features.plot = "CC2", group.by = "dataset", do.return = TRUE, point.size.use = 0)
p3 <- VlnPlot(object = hypo, features.plot = "CC3", group.by = "dataset", do.return = TRUE, point.size.use = 0)


pdf("Figures/CCA_plots.pdf")
plot_grid(p, p1, p2, p3)
dev.off()

PrintDim(object = hypo, reduction.type = "cca", dims.print = 1:2, genes.print = 10)

# Find significant CCAs for reduction using heatmaps

pdf("Figures/CCA_1-9.pdf")
DimHeatmap(object = hypo, 
           reduction.type = "cca", 
           cells.use = 500, 
           dim.use = 1:9, 
           do.balanced = TRUE)
dev.off()

pdf("Figures/CCA_10-18.pdf")
DimHeatmap(object = hypo, 
           reduction.type = "cca", 
           cells.use = 500, 
           dim.use = 10:18, 
           do.balanced = TRUE)
dev.off()

pdf("Figures/CCA_19-27.pdf")
DimHeatmap(object = hypo, 
           reduction.type = "cca", 
           cells.use = 500, 
           dim.use = 19:27, 
           do.balanced = TRUE)
dev.off()

pdf("Figures/CCA_28-36.pdf")
DimHeatmap(object = hypo, 
           reduction.type = "cca", 
           cells.use = 500, 
           dim.use = 28:36, 
           do.balanced = TRUE)
dev.off()

pdf("Figures/CCA_37-48.pdf")
DimHeatmap(object = hypo, 
           reduction.type = "cca", 
           cells.use = 500, 
           dim.use = 37:48, 
           do.balanced = TRUE)
dev.off()

# Remove cells who's differences aren't explained by CCA

hypo <- CalcVarExpRatio(object = hypo, reduction.type = "pca", grouping.var = "dataset", 
                        dims.use = 1:50)

# We discard cells where the variance explained by CCA is <2-fold (ratio <
# 0.5) compared to PCA
hypo.all.save <- hypo
hypo <- SubsetData(object = hypo, subset.name = "var.ratio.pca", accept.low = 0.5)

hypo.discard <- SubsetData(object = hypo.all.save, subset.name = "var.ratio.pca", 
                           accept.high = 0.5)
median(x = hypo@meta.data[, "nGene"])

median(x = hypo.discard@meta.data[, "nGene"])

# This removes ~750 cells from analysis, not bad
dim(hypo.discard@meta.data)
# [1] 750   7
dim(hypo@meta.data)
# [1] 17225     7
dim(hypo.all.save@meta.data)
# [1] 17975     5

rm(hypo.all.save, hypo.discard) # remove extra objects

# Next, align subsace using calculated CCAs

hypo <- AlignSubspace(object = hypo, reduction.type = "cca", grouping.var = "dataset", 
                      dims.align = 1:50)

# Visualize the aligned CCA and perform integrated analysis

pa <- DimPlot(object = hypo, reduction.use = "cca", group.by = "dataset", pt.size = 0.5, 
             do.return = TRUE)

pa1 <- VlnPlot(object = hypo, features.plot = "ACC1", group.by = "dataset", 
              do.return = TRUE, point.size.use = 0)
pa2 <- VlnPlot(object = hypo, features.plot = "ACC2", group.by = "dataset", 
              do.return = TRUE, point.size.use = 0)
pa3 <- VlnPlot(object = hypo, features.plot = "ACC3", group.by = "dataset", 
               do.return = TRUE, point.size.use = 0)

pdf("Figures/CCA_plots_aligned.pdf")
plot_grid(pa, pa1, pa2, pa3)
dev.off()


# Run analysis on all cells

hypo <- RunTSNE(object = hypo, reduction.use = "cca.aligned", dims.use = 1:50, 
                do.fast = TRUE)
hypo <- FindClusters(object = hypo, reduction.type = "cca.aligned", dims.use = 1:50, 
                     save.SNN = TRUE)
p1 <- TSNEPlot(object = hypo, group.by = "dataset", do.return = TRUE, pt.size = 0.5)
p2 <- TSNEPlot(object = hypo, do.return = TRUE, pt.size = 0.5)
p3 <- TSNEPlot(object = hypo, group.by = "res.1.2", do.return = TRUE, pt.size = 0.5)
p4 <- TSNEPlot(object = hypo, group.by = "res.0.6", do.return = TRUE, pt.size = 0.5)

pdf("Figures/CCA_plots_cluster_CCA.pdf", width = 16, height = 14)
plot_grid(p1, p2, p3, p4)
dev.off()

pdf("Figures/hypo_cluster_feature.plot_1.pdf")
FeaturePlot(object = hypo, features.plot = c("snap25a", "snap25b", "syt1a"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = .5)
dev.off()

# Calculate PCA and TSNE plot for PCA

hypo <- RunPCA(object = hypo, pc.genes = var.union, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 100)

hypo <- RunTSNE(object = hypo, dims.use = 1:50, check_duplicates = FALSE, do.fast = TRUE)


p1 <- TSNEPlot(object = hypo, group.by = "dataset", do.return = TRUE, pt.size = 0.5)
p2 <- TSNEPlot(object = hypo, do.return = TRUE, pt.size = 0.5)
p3 <- TSNEPlot(object = hypo, group.by = "res.1.2", do.return = TRUE, pt.size = 0.5)
p4 <- TSNEPlot(object = hypo, group.by = "res.0.6", do.return = TRUE, pt.size = 0.5)

pdf("Figures/CCA_plots_cluster_PCA.pdf", width = 16, height = 14)
plot_grid(p1, p2, p3, p4)
dev.off()

pdf("Figures/hypo_cluster_feature.plot_1_pca.pdf")
FeaturePlot(object = hypo, features.plot = c("snap25a", "snap25b", "syt1a"), min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = .5)
dev.off()

# Tables of mapping of clusters + figures to compare to random Forest

table(hypo@meta.data$dataset, hypo@meta.data$res.0.8)

#           0    1   10   11   12   13   14   15   16   17   18   19    2   20   21   22   23   24   25   26   27   28   29    3    4    5    6    7    8    9
# Raj    1033  228   41   65   30   24   79   62  125   22   18   27  279   10  104   18   20   14    4   18    2    7   20  123  127  153  157  186  147   37
# Shafer 3133 1122  421  353  345  345  287  278  182  281  249  224  900  223  111  155  147  147   92   55   69   59   42  962  912  699  655  598  530  469

x <- as.matrix(table(hypo@meta.data$res.0.6, hypo@meta.data$res.0.8))
row.names(x) <- paste("Jv", row.names(x), sep = "")
colnames(x) <- as.character(colnames(x))

y <- as.matrix(table(hypo@meta.data$res.1.2, hypo@meta.data$res.0.8))
row.names(y) <- paste("Ad", row.names(y), sep = "")
colnames(y) <- as.character(colnames(y))

z <- rbind(x, y)

pdf("Figures/Confusion_matrix_seurat_clustering.pdf", height = 16, width = 14)
plotConfusionMatrix(z,
                    x.order = as.character(c(0:29)),
                    y.order = rev(JvAd),
                    row.scale=TRUE,
                    col.scale=FALSE,
                    max.size = 12, 
                    xlab.use="CCA_Clusters", 
                    ylab.use="Combined_PCA_Clusters") 
dev.off()

pdf("Figures/Confusion_matrix_seurat_clustering_adult.pdf", height = 12, width = 14)
plotConfusionMatrix(table(hypo@meta.data$res.1.2, hypo@meta.data$res.0.8), 
                    x.order = as.character(c(0:29)), 
                    y.order = as.character(c(36:0)), 
                    row.scale=TRUE, 
                    max.size = 12, 
                    xlab.use="CCA_Clusters", 
                    ylab.use="Shafer_PCA_Clusters") 
dev.off()

pdf("Figures/Confusion_matrix_seurat_clustering_juvenile.pdf", height = 6, width = 14)
plotConfusionMatrix(table(hypo@meta.data$res.0.6, hypo@meta.data$res.0.8), 
                    x.order = as.character(c(0:29)), 
                    y.order = as.character(c(11:0)), 
                    row.scale=TRUE, 
                    max.size = 12, 
                    xlab.use="CCA_Clusters", 
                    ylab.use="Raj_PCA_Clusters") 
dev.off()

# Generate feature plots for expression of original cluster markers (shafer PCA clusters)

for (i in 0:36) {
  cluster <- hypo.markers %>% filter(cluster == i) %>% top_n(10, avg_logFC)
  pdf(paste("Figures/cluster_markers/hypo.cluster.markers.", i, ".pdf", sep = ""))
  FeaturePlot(object = hypo, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.1)
  dev.off()
}


hypo_combined_cca <- hypo

save(hypo_combined_cca, file = "hypo_combined_cca.Robj")
