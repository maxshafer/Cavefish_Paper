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

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

## Load objects
hypo <- readRDS("DanRer_65k.rds")

DefaultAssay(hypo) <- "RNA"

# Make Dendrograms as a new list

dendrograms <- list()

Idents(hypo) <- "Subtype"
hypo <- BuildClusterTree(hypo)
dendrograms[[1]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "SubclusterType"
hypo <- BuildClusterTree(hypo)
dendrograms[[2]] <- Tool(hypo, slot = "BuildClusterTree")

names(dendrograms) <- c("Subtype", "SubclusterType")

saveRDS(dendrograms, file = "Dan_dendrograms.rds")
