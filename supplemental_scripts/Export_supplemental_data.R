library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(reshape2)

## Load objects for exporting meta.data

hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")
hypo.integrated <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Hypo_integrated_127k_1500VFs_100Dims_v3.rds")

hypo.zeb.meta <- hypo.zeb@meta.data
colnames(hypo.zeb.meta)[7:9] <- c("Cluster", "Subcluster", "Subcluster_number")
hypo.ast.meta <- hypo.ast@meta.data
colnames(hypo.ast.meta)[10:16] <- c("Cluster", "Subcluster", "Subcluster_number", "species_Cluster", "species_Subcluster", "morph_Cluster", "morph_Subcluster")
hypo.integrated.meta <- hypo.integrated@meta.data
colnames(hypo.integrated.meta)[c(7:9,13:16,21,23,24,26,27)] <- c("Cluster", "Subcluster", "Subcluster_number", "species_Cluster", "species_Subcluster", "morph_Cluster", "morph_Subcluster", "integrated_Cluster", "integrated_Subcluster", "integrated_Subcluster_number", "integrated_Cluster_species", "integrated_Subcluster_species")

# Export to csv
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/object_metadata/")

write.csv(hypo.zeb.meta, file = "Danio_rerio_meta_data.csv")
write.csv(hypo.ast.meta, file = "Astyanax_mexicanus_meta_data.csv")
write.csv(hypo.integrated.meta, file = "Integrated_meta_data.csv")


## Load normed data for saving as csv

norm.cluster <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Normed_expression_data.rds")

names(norm.cluster) <- c("Clusters", "Subclusters", "integrated_Clusters", "integrated_Subclusters")

for (i in 1:length(norm.cluster)) {
  names(norm.cluster[[i]]) <- c("Danio_rerio", "Astyanax_mexicanus", "Astyanax_mexicanus_surface", "Astyanax_mexcianus_cave")
}

# Export to csv

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/pseudobulk_expression/")

for (i in 1:length(norm.cluster)) {
  for (j in 1:length(norm.cluster[[i]])) {
    write.csv(norm.cluster[[i]][[j]], file = paste(names(norm.cluster)[i], "-", names(norm.cluster[[i]])[j], ".csv", sep = ""))
  }
}

## Load TF_RF Data for integrated and ast
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC")
gene.lists <- c("nps", "nts", "synaptic", "ion")

gene.lists.int <- list()
gene.lists.int$Link_lists_Danio_rerio <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
gene.lists.int$Link_lists_Astyanax_mexicanus <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
gene.lists.int$tfModules_Danio_rerio <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
gene.lists.int$tfModules_Astyanax_mexicanus <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))

for (i in 1:length(gene.lists.int)) {
  names(gene.lists.int[[i]]) <- c("Neuropeptides", "Neurotransmitters", "Synaptic_genes", "Ion_channels")
}

# Export to csv

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/SCENIC_results/")

for (i in 1:length(gene.lists.int)) {
  for (j in 1:length(gene.lists.int[[i]])) {
    write.csv(gene.lists.int[[i]][[j]], file = paste(names(gene.lists.int)[i], "-", names(gene.lists.int[[i]])[j], ".csv", sep = ""))
  }
}

# Load GENIE3 files from all the comps, for ast and zeb

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC")
gene.lists <- c("nps")

gene.lists.ast <- list()
gene.lists.ast$Link_lists_surface <- lapply(gene.lists, function(x) readRDS(paste("surface/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
gene.lists.ast$Link_lists_cave <- lapply(gene.lists, function(x) readRDS(paste("cave/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
gene.lists.ast$tfModules_surface <- lapply(gene.lists, function(x) readRDS(paste("surface/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
gene.lists.ast$tfModules_cave <- lapply(gene.lists, function(x) readRDS(paste("cave/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))

for (i in 1:length(gene.lists.ast)) {
  names(gene.lists.ast[[i]]) <- c("Neuropeptides")
}

# Export to csv

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/SCENIC_results/")

for (i in 1:length(gene.lists.int)) {
  for (j in 1:length(gene.lists.int[[i]])) {
    write.csv(gene.lists.int[[i]][[j]], file = paste(names(gene.lists.int)[i], "-", names(gene.lists.int[[i]])[j], ".csv", sep = ""))
  }
}

## Load weir fst genes

weir.genes <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/astyanax_variants/weir.genes.rds")

# Export to csv

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Supplemental_data/weir_genes/")

for (i in 1:length(weir.genes)) {
  write.csv(weir.genes[[i]], file = paste(names(weir.genes)[i], ".csv", sep = ""))
}

