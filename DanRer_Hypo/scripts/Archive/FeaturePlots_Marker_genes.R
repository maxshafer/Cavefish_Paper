library(Seurat)
library(Matrix)
library(dplyr)

setwd("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo")

load("Shafer_Hypo_66k.Robj")

idents <- levels(hypo.zeb@ident)
hypo.zeb <- SetAllIdent(hypo.zeb, id = "res.0.6")
clusters <- levels(hypo.zeb@ident)

# Load and extract genes for plotting
# Make a list, named for each cluster

test <- read.csv("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/CSV/full_dataset/Cluster_annotations_Shafer_Hypo.res.0.6.csv", header = TRUE, row.names = "cluster")

genes <- strsplit(as.character(test$Marker.genes), "/")

names(genes) <- paste("Cluster", clusters, idents, sep = "_")

genes <- genes[3:length(genes)]

# Lapply a FeaturePlot png save
# x = element of genes list

saveFeature <- function(features = features, object = object, name = name) {
	png(filename = paste(name, "marker_genes.png", sep = "_"), units = "in", res = 250, height = 2, width = 4)
	FeaturePlot(object = object, features.plot = unlist(features), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, no.axes = T)
	dev.off()
}

setwd("~/Documents/Schier_Lab/R_Projects/Shafer_Hypo/Vassilis")

lapply(seq_along(genes), function(x) saveFeature(features = genes[[x]], object = hypo.zeb, name = names(genes)[[x]]))
