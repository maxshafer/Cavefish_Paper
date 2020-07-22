# library(networkD3)
library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?
Idents(hypo.integrated) <- "species.2"
hypo.ast <- subset(hypo.integrated, idents = "astyanax")

prop.table <- table(hypo.ast@meta.data$SubclusterType, hypo.ast@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$surface_specific <- ifelse(prop.table$astyanax_surface > .9, "yes", "no")
prop.table$cave_specific <- ifelse(prop.table$astyanax_surface < .1, "yes", "no")

surface.names <- row.names(prop.table[prop.table$surface_specific == "yes",])
cave.names <- row.names(prop.table[prop.table$cave_specific == "yes",])

# Choose subset of data
Idents(hypo.integrated) <- "integrated_Subtype"

hypo <- subset(hypo.integrated, idents = "GABA_2")

Idents(hypo.integrated) <- "integrated_SubclusterType"

# Find the integrated identities that a species specific subcluster maps too, for mapping the connections
test <- hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax" & hypo.integrated@meta.data$SubclusterType %in% cave.names[[1]], "integrated_SubclusterType"]
test <- table(test)
test <- test[test > sum(test)/10]

hypo <- subset(hypo.integrated, idents = names(test))

# Set up data

table1 <- table(hypo@meta.data[hypo@meta.data$species.2 == "astyanax", "Subtype"], hypo@meta.data[hypo@meta.data$species.2 == "astyanax","integrated_Subtype"])
table2 <- table(hypo@meta.data[hypo@meta.data$species.2 == "zebrafish", "integrated_Subtype"], hypo@meta.data[hypo@meta.data$species.2 == "zebrafish","Subtype"])

table1 <- table(hypo@meta.data[hypo@meta.data$species.2 == "astyanax", "SubclusterType"], hypo@meta.data[hypo@meta.data$species.2 == "astyanax","integrated_SubclusterType"])
table2 <- table(hypo@meta.data[hypo@meta.data$species.2 == "zebrafish", "integrated_SubclusterType"], hypo@meta.data[hypo@meta.data$species.2 == "zebrafish","SubclusterType"])

table1 <- melt(table1)
table1$Var.1 <- paste("ast", table1$Var.1, sep = "_")
table1$Var.2 <- paste("int", table1$Var.2, sep = "_")

table2 <- melt(table2)
table2$Var.1 <- paste("int", table2$Var.1, sep = "_")
table2$Var.2 <- paste("zeb", table2$Var.2, sep = "_")

nodes <- c(as.character(unique(table1$Var.2)), as.character(unique(table1$Var.1)), as.character(unique(table2$Var.2)))
nodes <- as.data.frame(nodes)


table <- rbind(table1, table2)

table.ready <- table[table$value > 10,]

## Plot, as html file (can save using RStudio)

nodes2 <- nodes[nodes$nodes %in% unique(c(table.ready$Var.1, table.ready$Var.2)),]
nodes2 <- as.data.frame(nodes2)

table.ready$Var1 <- as.numeric(plyr::mapvalues(x = table.ready$Var.1, from = nodes2$nodes2, to = row.names(nodes2)))-1
table.ready$Var2 <- as.numeric(plyr::mapvalues(x = table.ready$Var.2, from = nodes2$nodes2, to = row.names(nodes2)))-1

# Plot sankeyNetwork
sankeyNetwork(Links = table.ready, Nodes = nodes2, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes2", units = "cells", fontSize = 8, nodeWidth = 40, nodePadding = 10, height = 200, width = 300, iterations = 10000)

sankeyNetwork(Links = table.ready, Nodes = nodes2, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes2", units = "cells", fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 800, width = 800, iterations = 10000)

sum(table.ready$value)
sum(table.ready$value[grep("zeb", table.ready$Var.2)])

