library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(reshape2)

# Load hypo and make species specific names (for subsetting ast gene lists)

hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")

Idents(hypo) <- "SubclusterType"

# Find morph specific clusters
prop.table <- table(hypo@meta.data$SubclusterType, hypo@meta.data$species)
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo))
index <- c(surface.names, cave.names)
index2 <- subclusters[!(subclusters %in% c(surface.names, cave.names))]


# Load gene lists for clusters and subclusters

gene.lists.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/marker_gene_lists.rds")

gene.lists.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/marker_gene_lists.rds")

gene.lists <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/drift_gene_lists.rds")


# Adjust gene.lists.ast

gene.lists.ast$surface.specific.markers.sub <- gene.lists.ast$surface.markers.sub[index[index %in% names(gene.lists.ast$surface.markers.sub)]]
gene.lists.ast$cave.specific.markers.sub <- gene.lists.ast$cave.specific.markers.sub[index[index %in% names(gene.lists.ast$cave.specific.markers.sub)]]

gene.lists.ast$conserved.markers.sub <- gene.lists.ast$conserved.markers.sub[index2[index2 %in% names(gene.lists.ast$conserved.markers.sub)]]
gene.lists.ast$surface.markers.sub <- gene.lists.ast$surface.markers.sub[index2[index2 %in% names(gene.lists.ast$surface.markers.sub)]]
gene.lists.ast$cave.markers.sub <- gene.lists.ast$cave.markers.sub[index2[index2 %in% names(gene.lists.ast$cave.markers.sub)]]


# Add cluster/subcluster column to each data.frame

for (i in 1:length(gene.lists.zeb)) {
  for (j in 1:length(gene.lists.zeb[[i]])) {
    gene.lists.zeb[[i]][[j]]$cluster <- names(gene.lists.zeb[[i]])[j]
  }
}

for (i in 1:length(gene.lists.ast)) {
  for (j in 1:length(gene.lists.ast[[i]])) {
    gene.lists.ast[[i]][[j]]$cluster <- names(gene.lists.ast[[i]])[j]
  }
}

for (i in 1:length(gene.lists)) {
  for (j in 1:length(gene.lists[[i]])) {
    gene.lists[[i]][[j]]$cluster <- names(gene.lists[[i]])[j]
  }
}

# Concatenate data.frames

gene.lists.zeb.2 <- lapply(seq_along(gene.lists.zeb), function(x) Reduce(rbind, gene.lists.zeb[[x]]))
names(gene.lists.zeb.2) <- c("zebrafish.markers", "zebrafish.markers.sub")

gene.lists.ast.2 <- lapply(seq_along(gene.lists.ast), function(x) Reduce(rbind, gene.lists.ast[[x]])) # df 4 doesn't work b/c the morph specific clusters have different dimensions
names(gene.lists.ast.2) <- names(gene.lists.ast)

gene.lists.2 <- lapply(seq_along(gene.lists), function(x) Reduce(rbind, gene.lists[[x]]))
names(gene.lists.2) <- names(gene.lists)

# Export to csv

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/marker_gene_lists/")

for (i in 1:length(gene.lists.zeb.2)) {
  write.csv(gene.lists.zeb.2[[i]], file = paste("Drerio_", names(gene.lists.zeb.2)[i], ".csv", sep = ""))
}

for (i in 1:length(gene.lists.ast.2)) {
  write.csv(gene.lists.ast.2[[i]], file = paste("Amexicanus_", names(gene.lists.ast.2)[i], ".csv", sep = ""))
}

for (i in 1:length(gene.lists.2)) {
  write.csv(gene.lists.2[[i]], file = paste("Integrated_", names(gene.lists.2)[i], ".csv", sep = ""))
}

## Saves all to directory, then can zip and upload as one supplemental data file







