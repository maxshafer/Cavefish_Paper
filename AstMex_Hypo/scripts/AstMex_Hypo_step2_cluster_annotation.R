library(Seurat)
library(Matrix)
library(plyr)
library(stringr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo")

## Load objects

test <- read.csv("CSV/full_dataset/Cluster_annotations_AstMex_Hypo.res.0.6.csv", header = TRUE, row.names = "cluster")
load("Archive/AstMex_66k.Robj")
hypo <- hypo.ast
rm(hypo.ast)

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

## Assign subcluster ids by paste
## subcluster.NEW doesn't exist anymore?? whered it go? res.0.4?

hypo@meta.data$SubclusterType <- paste(hypo@meta.data$Subtype, sub(".*_", "", hypo@meta.data$subcluster.NEW), sep = "_")

hypo.ast <- hypo

save(hypo.ast, file = "AstMex_65k.Robj")

## Remove contaminating clusters and save as different object
## One whole subtype, plus one subcluster

hypo <- SetAllIdent(hypo, id = "Subtype")

# Remove all 3 CONT, or keep everything else

hypo.ast <- SubsetData(hypo, ident.use = levels(hypo@ident)[!grepl("CONT", levels(hypo@ident))])

hypo.ast <- SetAllIdent(hypo.ast, id = "SubclusterType")
hypo.ast <- SubsetData(hypo.ast, ident.use = levels(hypo.ast@ident)[!grepl("GABA_5_1", levels(hypo.ast@ident))])

# Rename immune cell types

hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_1"] <- "Bcells"
hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_2"] <- "Mast_cells"
hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_3"] <- "Thrombocytes"
hypo.ast@meta.data$Subtype[hypo.ast@meta.data$Subtype == "Leucocytes_4"] <- "Neutrophils"

# Set factor levels for Subtypes

hypo.ast@meta.data$Subtype <- factor(hypo.ast@meta.data$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated/Ventricular", "Ependymal", "Progenitors_1", "Progenitors_2", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes_1", "Oligodendrocytes_2", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "GABA_6", "GABA_7", "GABA_Prdx1_1", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Galanin", "Otpa/b_1", "Otpa/b_2", "Otpa/b_3", "Lymphatic", "Tcells", "Bcells", "Mast_cells", "Thrombocytes", "Neutrophils", "Macrophages", "Microglia"))

## Assign subcluster ids by paste

hypo.ast@meta.data$SubclusterType <- paste(hypo.ast@meta.data$Subtype, sub(".*_", "", hypo.ast@meta.data$res.0.4), sep = "_")

## Set factor levels for SubclusterTypes based on Subtype factors levels

hypo.ast@meta.data$SubclusterType <- factor(hypo.ast@meta.data$SubclusterType, levels = unique(hypo.ast@meta.data$SubclusterType)[unlist(sapply(levels(hypo.ast@meta.data$Subtype), function(x) grep(x, unique(hypo.ast@meta.data$SubclusterType))))])

# Save as new, reduced file

save(hypo.ast, file = "AstMex_64k.Robj")



