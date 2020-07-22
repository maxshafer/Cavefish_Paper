library(networkD3)
library(alluvial)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

# Choose subset of data
Idents(hypo.integrated) <- "integrated_Subtype"

hypo <- subset(hypo.integrated, idents = "Glut_1")

# Set up data

table1 <- table(hypo@meta.data[hypo@meta.data$species.2 == "astyanax", "Subtype"], hypo@meta.data[hypo@meta.data$species.2 == "astyanax","integrated_Subtype"])
table2 <- table(hypo@meta.datta[hypo@meta.data$species.2 == "zebrafish", "integrated_Subtype"], hypo@meta.data[hypo@meta.data$species.2 == "zebrafish","Subtype"])


table1 <- melt(table1)
table2 <- melt(table2)

nodes <- c(as.character(unique(table1[,2])), as.character(unique(table1[,1])), as.character(unique(table2[,2])))
nodes <- as.data.frame(nodes)

# Make sure that integrated types appear as sources for one, and targets for the other species
table1$Var.2 <- as.numeric(table1$Var.2)-1
table1$Var.1 <- as.numeric(table1$Var.1)+max(table1$Var.2)
table2$Var.1 <- as.numeric(table2$Var.1)-1
table2$Var.2 <- as.numeric(table2$Var.2)+max(table1$Var.1)

table <- rbind(table1, table2)

table.ready <- table[table$value > 100,]

## Plot, as html file (can save using RStudio)

sankeyNetwork(Links = table.ready, Nodes = nodes, Source = "Var.1", Target = "Var.2", Value = "value", NodeID = "nodes", units = "cells", fontSize = 12, nodeWidth = 80, nodePadding = 5, height = 750, width = 750, iterations = 10000)
