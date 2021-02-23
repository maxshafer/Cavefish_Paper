library(networkD3)
library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

# Set up data

table1 <- table(hypo@meta.data[hypo@meta.data$species.2 == "astyanax", "Cluster"], hypo@meta.data[hypo@meta.data$species.2 == "astyanax","integrated_Cluster"])
table2 <- table(hypo@meta.data[hypo@meta.data$species.2 == "zebrafish", "integrated_Cluster"], hypo@meta.data[hypo@meta.data$species.2 == "zebrafish","Cluster"])

table1 <- melt(table1)
table1$Var.1 <- paste("ast", table1$Var.1, sep = "_")
table1$Var.2 <- paste("int", table1$Var.2, sep = "_")

table2 <- melt(table2)
table2$Var.1 <- paste("int", table2$Var.1, sep = "_")
table2$Var.2 <- paste("zeb", table2$Var.2, sep = "_")

nodes <- c(as.character(unique(table1$Var.2)), as.character(unique(table1$Var.1)), as.character(unique(table2$Var.2)))
nodes <- as.data.frame(nodes)


table <- rbind(table1, table2)

table.ready <- table[table$value > 100,]

## Plot, as html file (can save using RStudio)

nodes2 <- nodes[nodes$nodes %in% unique(c(table.ready$Var.1, table.ready$Var.2)),]
nodes2 <- as.data.frame(nodes2)
nodes2$Var <- as.numeric(row.names(nodes2))-1

table.ready$Var1 <- nodes2$Var[match(table.ready$Var.1, nodes2$nodes2)]
table.ready$Var2 <- nodes2$Var[match(table.ready$Var.2, nodes2$nodes2)]


# Plot sankeyNetwork
# "Show in new window" to show in browswer, then save as PDF for import
# Colours were changed in affinity

sankeyNetwork(Links = table.ready, Nodes = nodes2, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes2", units = "cells", fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 800, width = 400, iterations = 10000, sinksRight = F)


