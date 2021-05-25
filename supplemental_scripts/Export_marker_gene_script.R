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

a = 1.5
b = 2
f = 0.1

gene.lists.pos <- readRDS(paste("drift_gene_lists_pos_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

## Export the number of conserved and species specific marker genes per cluster

SI.numbers.species <- lapply(seq_along(1:length(gene.lists.pos[[7]])), function(x) {
  Zebrafish_specific <- length(setdiff(rownames(gene.lists.pos[[8]][[x]]), rownames(gene.lists.pos[[7]][[x]])))
  Mexicantetra_specific <- length(setdiff(rownames(gene.lists.pos[[9]][[x]]), rownames(gene.lists.pos[[7]][[x]])))
  Conserved <- length(rownames(gene.lists.pos[[7]][[x]]))
  return(data.frame(Conserved, Mexicantetra_specific, Zebrafish_specific)) } )

SI.numbers.species <- Reduce(rbind, SI.numbers.species)
SI.numbers.species$cell_type <- names(gene.lists.pos[[7]])

SI.numbers.morphs <- lapply(seq_along(1:length(gene.lists.pos[[12]])), function(x) {
  Surface_specific <- length(setdiff(rownames(gene.lists.pos[[13]][[x]]), rownames(gene.lists.pos[[12]][[x]])))
  Cave_specific <- length(setdiff(rownames(gene.lists.pos[[14]][[x]]), rownames(gene.lists.pos[[12]][[x]])))
  Conserved <- length(rownames(gene.lists.pos[[12]][[x]]))
  return(data.frame(Conserved, Cave_specific, Surface_specific)) } )

SI.numbers.morphs <- Reduce(rbind, SI.numbers.morphs)
SI.numbers.morphs$cell_type <- names(gene.lists.pos[[12]])

# Export to csv

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/SI_numbers_lists/")

write.csv(SI.numbers.species, file = "SI_numbers_species.csv")
write.csv(SI.numbers.morphs, file = "SI_numbers_morphs.csv")


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







