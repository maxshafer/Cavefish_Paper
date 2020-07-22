library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dendextend)
library(ape)
library(circlize)
library(ggplot2)
library(scales)
library(phytools)
library(ggtree)
library(ggplotify)

# Load seurat object (astmex)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k.rds")
DefaultAssay(hypo) <- "RNA"


# Find morph specific clusters
prop.table <- table(hypo@meta.data$SubclusterType, hypo@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$surface_specific <- ifelse(prop.table$astyanax_surface > .9, "yes", "no")
prop.table$cave_specific <- ifelse(prop.table$astyanax_cave > .9, "yes", "no")

surface.names <- row.names(prop.table[prop.table$surface_specific == "yes",])
cave.names <- row.names(prop.table[prop.table$cave_specific == "yes",])

# Find cave specific clusters
prop.table <- table(hypo@meta.data$SubclusterType, hypo@meta.data$morph)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$surface_specific <- ifelse(prop.table$Choy_surface > .9, "yes", "no")
prop.table$Molino_specific <- ifelse(prop.table$Molino_cave > .9, "yes", "no")
prop.table$Pachon_specific <- ifelse(prop.table$Pachon_cave > .9, "yes", "no")
prop.table$Tinaja_specific <- ifelse(prop.table$Tinaja_cave > .9, "yes", "no")

surface.names <- row.names(prop.table[prop.table$surface_specific == "yes",])
cave.names <- row.names(prop.table[prop.table$cave_specific == "yes",])

# Make index
subclusters <- levels(Idents(hypo))
index <- subclusters[!(subclusters %in% c(surface.names, cave.names))]

# Make Dendrograms as a new list

dendrograms <- list()

Idents(hypo) <- "Subtype"
hypo <- BuildClusterTree(hypo)
dendrograms[[1]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "SubclusterType"
hypo <- BuildClusterTree(hypo)
dendrograms[[2]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "species_Subtype"
hypo <- BuildClusterTree(hypo)
dendrograms[[3]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "SubclusterType"
hypo.1 <- subset(hypo, idents = index)
Idents(hypo.1) <- "species_SubclusterType"
hypo.1 <- BuildClusterTree(hypo.1)
dendrograms[[4]] <- Tool(hypo.1, slot = "BuildClusterTree")

Idents(hypo) <- "morph_Subtype"
hypo <- BuildClusterTree(hypo)
dendrograms[[5]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "morph_SubclusterType"
hypo <- BuildClusterTree(hypo)
dendrograms[[6]] <- Tool(hypo, slot = "BuildClusterTree")

saveRDS(dendrograms, file = "Ast_dendrograms.rds")
