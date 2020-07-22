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

names (data.vec) <- c("hypo1","hypo2","hypo3","hypo4",
					  "hypo5","hypo6","hypo7","hypo8",
					  "hypo9","hypo10","hypo11","hypo12",
					  "hypo13","hypo14","hypo15","hypo16")

input.data <- Read10X(data.vec)

# Create Seurat object, normalize, scale and find variable genes

hypo <- CreateSeuratObject(input.data, 
							min.genes = 200, 
							project = "Shafer_Hypo") # 32191 genes across 66993 samples.

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

hypo@meta.data$sex2[hypo@meta.data$sex == "male"] <- "male"
hypo@meta.data$sex2[hypo@meta.data$sex == "male2"] <- "male"
hypo@meta.data$sex2[hypo@meta.data$sex == "female1"] <- "female"
hypo@meta.data$sex2[hypo@meta.data$sex == "female2"] <- "female"


# QC figures

pdf("Figures/VlnPlot_nGene_nUMI_orig.ident.pdf")
VlnPlot(hypo, 
		features.plot = c("nGene", "nUMI"), 
		group.by = "orig.ident", 
		nCol = 2, 
		point.size.use = 0.01)
dev.off()

pdf("Figures/VlnPlot_nGene_nUMI_sex.pdf")
VlnPlot(hypo, 
		features.plot = c("nGene", "nUMI"), 
		group.by = "sex", 
		nCol = 2, 
		point.size.use = 0.01)
dev.off()

pdf("Figures/GenePlot_nGene_nUMI.pdf")
GenePlot(hypo, 
		gene1 = "nUMI", 
		gene2 = "nGene")
dev.off()

# Filter out cells with nGene higher then 3k, nUMI higher then 15k?

hypo <- FilterCells(object = hypo, 
					subset.names = c("nGene"), 
					low.thresholds = c(200), 
					high.thresholds = c(2500)) # 32191 genes across 65947 samples.

# Normalize, scale and find variable genes

hypo <- NormalizeData(object = hypo, 
						normalization.method = "LogNormalize", 
						scale.factor = 10000)
						
pdf("Figures/VariableGenes.pdf")
hypo <- FindVariableGenes(object = hypo, 
						mean.function = ExpMean, 
						dispersion.function = LogVMR, 
						do.plot = TRUE)
dev.off()

hypo <- ScaleData(object = hypo, 
				genes.use = hypo@var.genes, 
				do.par = TRUE)
				

# Run PCA for clustering

hypo <- RunPCA(object = hypo, 
				pc.genes = hypo@var.genes, 
				do.print = FALSE, 
				pcs.compute = 150)

# Determine the PCAs to use for clustering - using all 20

pdf("Figures/PCElbowPlot.pdf", height = 14, width = 14)
PCElbowPlot(object = hypo, 
			num.pc = 150)
dev.off()

pdf("Figures/PCHeatmap_1-5_45-50.pdf", height = 14, width = 14)
PCHeatmap(object = hypo, 
		pc.use = c(1:5, 45:50), 
		cellsuse = 500, 
		do.balanced = TRUE, 
		label.columns = FALSE)
dev.off()

pdf("Figures/PCHeatmap_1-5_70-75.pdf", height = 14, width = 14)
PCHeatmap(object = hypo, 
		pc.use = c(1:5, 70:75), 
		cellsuse = 500, 
		do.balanced = TRUE, 
		label.columns = FALSE)
dev.off()



# Find clusters! Run for multiple resolutions

hypo <- FindClusters(object = hypo, 
					reduction.type = "pca", 
					dims.use = 1:75, 
					resolution = 0.3, 
					save.SNN = TRUE, 
					n.start = 10, 
					nn.eps = 0.5, 
					print.output = FALSE, 
					force.recalc = TRUE)
					
hypo <- FindClusters(object = hypo, 
					resolution = 0.6)
hypo <- FindClusters(object = hypo, 
					resolution = 1.2)
hypo <- FindClusters(object = hypo, 
					resolution = 3.0)
hypo <- FindClusters(object = hypo, 
					resolution = 5.0)
hypo <- FindClusters(object = hypo, 
					resolution = 10.0)


# Calculate tSNE for plotting

set.seed(1188)
hypo <- RunTSNE(object = hypo, 
				reduction.use = "pca", 
				dims.use = 1:75, 
				tsne.method = "FIt-SNE", 
				nthreads = 6, 
				reduction.name = "FItSNE", 
				reduction.key = "FItSNE_", 
				fast_tsne_path = "~/Documents/Programs/FIt-SNE-master/bin/fast_tsne", 
				max_iter = 2000)

cols <- c("firebrick1", "indianred1", "royalblue1", "skyblue1")

cols2 <- c("royalblue1", "skyblue2", "skyblue3", "skyblue4", 
			"firebrick1", "firebrick2", "firebrick3", "firebrick4", 
			"royalblue2", "royalblue3", "royalblue4", "indianred1", 
			"indianred2", "indianred3", "indianred"," skyblue1")

cluster <- DimPlot(hypo, 
				reduction.use = "FItSNE", 
				do.label = TRUE, 
				no.legend = TRUE, 
				no.axes = TRUE, 
				do.return = TRUE, 
				vector.friendly = TRUE, pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))

sample <- DimPlot(hypo, 
				reduction.use = "FItSNE", 
				no.legend = FALSE, 
				no.axes = TRUE, 
				group.by = "orig.ident", 
				do.return = TRUE, 
				vector.friendly = TRUE, 
				pt.size = 0.5, 
				cols.use = cols2) + ggtitle("Sample") + theme(plot.title = element_text(hjust = 0.5))
				
sex <- DimPlot(hypo, 
				reduction.use = "FItSNE", 
				no.legend = TRUE, 
				no.axes = TRUE, 
				group.by = "sex", 
				do.return = TRUE, 
				vector.friendly = TRUE, 
				pt.size = 0.5, 
				cols.use = cols) + ggtitle("Sex") + theme(plot.title = element_text(hjust = 0.5))

plot_grid(cluster, sample, sex)


# Plot clusters for every resolution

pdf("Figures/hypo_cluster_plot_res.0.3.pdf", height = 7, width = 7)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID - res 0.3") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figures/hypo_cluster_plot_res.0.6.pdf", height = 7, width = 7)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "res.0.6", 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figures/hypo_cluster_plot_res.1.2.pdf", height = 7, width = 7)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "res.1.2", 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figures/hypo_cluster_plot_res.3.pdf", height = 7, width = 7)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "res.3", 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figures/hypo_cluster_plot_res.5.pdf", height = 7, width = 7)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "res.5", 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figures/hypo_cluster_plot_res.10.pdf", height = 7, width = 7)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "res.10", 
		do.return = TRUE, 
		vector.friendly = TRUE,
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf("Figures/hypo_cluster_plot_orig_ident.pdf", height = 7, width = 8)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "orig.ident", 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))dev.off()

pdf("Figures/hypo_cluster_plot_sex.pdf", height = 7, width = 8)
DimPlot(hypo, 
		reduction.use = "FItSNE", 
		do.label = TRUE, 
		no.legend = TRUE, 
		group.by = "sex", 
		do.return = TRUE, 
		vector.friendly = TRUE, 
		pt.size = 1) + ggtitle("Cluster ID") + theme(plot.title = element_text(hjust = 0.5))dev.off()


# VlnPlot(hypo, features.plot = c("thy1", "ecscr"), group.by = "sex", nCol = 1, point.size.use = 1)

# FeaturePlot(object = hypo, features.plot = c("nrn1a"), cols.use = c("grey85", "blue", "blue", "blue"), reduction.use = "FItSNE", min.cutoff = "q9", pt.size = .5, pch.use = 1)


freq_table <- prop.table(table(hypo@meta.data$res.0.3, hypo@meta.data$sex2), margin = 2)

identities <- levels(hypo@ident)
my_color_palette <- hue_pal()(length(identities))

pdf("Figures/hypo_cluster_sex_freq.pdf", height = 30, width = 30)
barplot(height = freq_table, legend.text = F, col = my_color_palette)
dev.off()


# Find markers

hypo.markers <- FindAllMarkers(object = hypo, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

# Generate Feature plot of top 2 markers each

features.markers <- hypo.markers %>% group_by(cluster) %>% top_n(1, avg_logFC)

pdf("Figures/hypo_cluster.markers.0.3.pdf", height = 60, width = 30)
FeaturePlot(object = hypo, features.plot = features.markers$gene, min.cutoff = "q9", cols.use = c("lightgrey", "blue"), pt.size = 0.5)
dev.off()

pdf("Figures/hypo_cluster.markers.0.3.vlnplot.pdf", height = 40, width = 20)
VlnPlot(object = hypo, features.plot = features.markers$gene, nCol = 1, size.x.use = 0, point.size.use = NA, size.title.use = 0)
dev.off()

# Generate lists + feature plots of top 10 markers for each cluster

for (i in 0:31) {
  name <- as.character(paste("cluster_", i, sep=""))
  assign(name, hypo.markers %>% filter(cluster == i) %>% top_n(20, avg_logFC))
  write.csv(name, file = paste("CSV/full_dataset/cluster_markers/cluster_", i, ".csv", sep=""))
}

# Generate feature plots for each cluster (top 20 genes)

for (i in 0:31) {
  cluster <- hypo.markers %>% filter(cluster == i) %>% top_n(9, avg_logFC)
  pdf(paste("Figures/full_dataset/cluster_markers/hypo.cluster.markers.0.3.", i, ".pdf", sep = ""), height = 15, width = 15)
  FeaturePlot(object = hypo, features.plot = cluster[,7], min.cutoff = "q9", cols.use = c("lightgrey", "blue"), nCol = 3, pt.size = 0.1)
  dev.off()
}


top10 <- hypo.markers.0.3 %>% group_by(cluster) %>% top_n(10, avg_logFC)
top20 <- hypo.markers.0.3 %>% group_by(cluster) %>% top_n(20, avg_logFC)

pdf("Figures/full_dataset/cluster.markers.1.2.heatmap.top10.pdf", height = 15, width = 100)
DoHeatmap(hypo, genes.use = top10$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()
pdf("Figures/full_dataset/cluster.markers.1.2.heatmap.top20.pdf", height = 30, width = 100)
DoHeatmap(hypo, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
dev.off()



# save seurat object and cluster markers

write.table(hypo.markers.1.2, file = "CSV/full_dataset/hypo.cluster.markers.1.2.txt", sep = "\t")

hypo.26k <- hypo
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



