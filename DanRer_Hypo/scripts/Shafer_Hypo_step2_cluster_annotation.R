library(Seurat)
library(Matrix)
library(plyr)
library(stringr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo")

## Load objects

test <- read.csv("CSV/full_dataset/Cluster_annotations_Shafer_Hypo.res.0.6.csv", header = TRUE, row.names = "cluster")
load("Shafer_Hypo_67k.Robj")
hypo <- hypo.zeb
rm(hypo.zeb)

hypo <- SetAllIdent(hypo, id = "res.0.6")
current.cluster.ids <- c(0:(length(levels(hypo@ident))-1))

new.cluster.ids <- as.vector(test$Subtype)
hypo@ident <- mapvalues(x = hypo@ident, from = current.cluster.ids, to = new.cluster.ids)
hypo <- StashIdent(hypo, save.name = "Subtype")

hypo <- SetAllIdent(hypo, id = "res.0.6")
new.cluster.ids <- as.vector(test$Type)
hypo@ident <- mapvalues(x = hypo@ident, from = current.cluster.ids, to = new.cluster.ids)
hypo <- StashIdent(hypo, save.name = "Type")

hypo <- SetAllIdent(hypo, id = "res.0.6")
new.cluster.ids <- as.vector(test$Type2)
hypo@ident <- mapvalues(x = hypo@ident, from = current.cluster.ids, to = new.cluster.ids)
hypo <- StashIdent(hypo, save.name = "Type2")

## Map new subcluster ids based on cluster merging

cluster_mapping <- read.csv("Cluster_mapping_combining.csv", stringsAsFactors = F)

hypo <- SetAllIdent(hypo, id = "subcluster.0.4")
current.cluster.ids <- as.numeric(levels(hypo@ident))

new.cluster.ids <- cluster_mapping$NEW
hypo@ident <- mapvalues(x = hypo@ident, from = current.cluster.ids, to = new.cluster.ids)
hypo <- StashIdent(hypo, save.name = "subcluster.NEW")

hypo@meta.data$SubclusterType <- paste(hypo@meta.data$Subtype, sub(".*_", "", hypo@meta.data$subcluster.NEW), sep = "_")

hypo.zeb <- hypo

save(hypo.zeb, file = "Shafer_Hypo_67k.Robj")

## Remove contaminating clusters and save as different object

hypo <- SetAllIdent(hypo, id = "Subtype")

# Remove all 3 CONT, or keep everything else

hypo.zeb <- SubsetData(hypo, ident.use = levels(hypo@ident)[!grepl("CONT", levels(hypo@ident))])

save(hypo.zeb, file = "Shafer_Hypo_66k.Robj")
