library(networkD3)
library(alluvial)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?

prop.table <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species.2)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$zebrafish_specific <- ifelse(prop.table$zebrafish > .90, "yes", "no")
prop.table$astyanax_specific <- ifelse(prop.table$astyanax > .90, "yes", "no")

zebrafish.names <- row.names(prop.table[prop.table$zebrafish_specific == "yes",])
astyanax.names <- row.names(prop.table[prop.table$astyanax_specific == "yes",])

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?
Idents(hypo.integrated) <- "species.2"
hypo.ast <- subset(hypo.integrated, idents = "astyanax")

prop.table <- table(hypo.ast@meta.data$SubclusterType, hypo.ast@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$surface_specific <- ifelse(prop.table$astyanax_surface > .9, "yes", "no")
prop.table$cave_specific <- ifelse(prop.table$astyanax_surface < .1, "yes", "no")

surface.names <- row.names(prop.table[prop.table$surface_specific == "yes",])
cave.names <- row.names(prop.table[prop.table$cave_specific == "yes",])

# Function to plot Sankey's for each

plotSankey <- function(object = hypo.integrated, species = "astyanax", cluster = "GABA_7_2", sum.div = 10, plot.cut = 10, remove.non = F, ...) {
  
  if (species == "astyanax") {
    species.2 <- "zebrafish"
  } else {
    species.2 <- "astyanax"
  }
  
  # Identify the integrated IDs that the target cluster maps to
  
  if (cluster %in% unique(object@meta.data$SubclusterType)) {
    res <- "SubclusterType"
    Idents(object) <- "integrated_SubclusterType"
    if (length(cluster) == 1) {
      table <- table(object@meta.data[object@meta.data$species.2 == species & object@meta.data$SubclusterType %in% cluster, "integrated_SubclusterType"])
      subset.names <- names(table[table > sum(table)/sum.div])
    }
    if (length(cluster) > 1) {
      subset.names <- list()
      for (i in 1:length(cluster)) {
        table <- table(object@meta.data[object@meta.data$species.2 == species & object@meta.data$SubclusterType %in% cluster[[i]], "integrated_SubclusterType"])
        subset.names[[i]] <- names(table[table > sum(table)/sum.div])
      }
      subset.names <- unlist(subset.names)
    }
  } 
  
  if (cluster %in% unique(object@meta.data$Subtype)) {
    res <- "Subtype"
    Idents(object) <- "integrated_Subtype"
    table <- table(object@meta.data[object@meta.data$species.2 == species & object@meta.data$Subtype %in% cluster, "integrated_Subtype"])
  }
  
  # Subset object to mapped integrated IDs
  #subset.names <- names(table[table > sum(table)/20]) # The problem: if using multiple idents, this should be normalized for the size of each
  print("Integrated identities used for plotting")
  print(subset.names)
  object <- subset(object, idents = subset.names)
  
  # Extract the count tables
  if (res == "SubclusterType") {
    if (!(species %in% object@meta.data$species.2)) {
      table2 <- melt(table(object@meta.data[object@meta.data$species.2 != species, "integrated_SubclusterType"], object@meta.data[object@meta.data$species.2 != species,"SubclusterType"]))
      table1 <- table2
      table1$value <- 0
      print(paste(cluster, "does not exist in", species, sep = " "))
    } else {
      if (!(species.2 %in% object@meta.data$species.2)) {
        table1 <- melt(table(object@meta.data[object@meta.data$species.2 == species, "SubclusterType"], object@meta.data[object@meta.data$species.2 == species,"integrated_SubclusterType"]))
        table2 <- table1
        table2$value <- 0
        print(paste(cluster, "does not exist in", species.2, sep = " "))
      } else {
        table1 <- melt(table(object@meta.data[object@meta.data$species.2 == species, "SubclusterType"], object@meta.data[object@meta.data$species.2 == species,"integrated_SubclusterType"]))
        table2 <- melt(table(object@meta.data[object@meta.data$species.2 != species, "integrated_SubclusterType"], object@meta.data[object@meta.data$species.2 != species,"SubclusterType"]))
      }
    }
  }
  if (res == "Subtype") {
    if (!(species %in% object@meta.data$species.2)) {
      table2 <- melt(table(object@meta.data[object@meta.data$species.2 != species, "integrated_Subtype"], object@meta.data[object@meta.data$species.2 != species,"Subtype"]))
      table1 <- table2
      table1$value <- 0
      print(paste(cluster, "does not exist in", species, sep = " "))
    } else {
      if (!(species.2 %in% object@meta.data$species.2)) {
        table1 <- melt(table(object@meta.data[object@meta.data$species.2 == species, "Subtype"], object@meta.data[object@meta.data$species.2 == species,"integrated_Subtype"]))
        table2 <- table1
        table2$value <- 0
        print(paste(cluster, "does not exist in", species.2, sep = " "))
      } else {
        table1 <- melt(table(object@meta.data[object@meta.data$species.2 == species, "Subtype"], object@meta.data[object@meta.data$species.2 == species,"integrated_Subtype"]))
        table2 <- melt(table(object@meta.data[object@meta.data$species.2 != species, "integrated_Subtype"], object@meta.data[object@meta.data$species.2 != species,"Subtype"]))
      }
    }
  }
  
  # Append species idents + combine lists + remove lowbounds
  table1$Var.1 <- paste(as.character(species), table1$Var.1, sep = "_")
  table1$Var.2 <- paste("int", table1$Var.2, sep = "_")
  table2$Var.1 <- paste("int", table2$Var.1, sep = "_")
  table2$Var.2 <- paste(as.character(species.2), table2$Var.2, sep = "_")
  
  table <- rbind(table1, table2)
  table <- table[table$value > plot.cut,]
  
  if (remove.non == TRUE) {
    cluster <- paste(species, cluster, sep = "_")
    table <- table[c(seq(1:length(table$Var.1))[table$Var.1 %in% cluster], grep("int", table$Var.1)),]
  }
  
  # Make nodes list
  nodes <- union(as.character(table$Var.2), as.character(table$Var.1))
  nodes <- as.data.frame(nodes)
  
  # Fix node levels
  table$Var1 <- as.numeric(plyr::mapvalues(x = table$Var.1, from = nodes$nodes, to = row.names(nodes)))-1
  table$Var2 <- as.numeric(plyr::mapvalues(x = table$Var.2, from = nodes$nodes, to = row.names(nodes)))-1
  
  return(sankeyNetwork(Links = table, Nodes = nodes, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes", units = "cells", iterations = 10000, sinksRight = F, ...))
}

## Make some figures

plotSankey(object = hypo.integrated, species = "astyanax", cluster = "GABA_7_2", sum.div = 10, plot.cut = 20, remove.non = F, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 600, width = 800)

plotSankey(object = hypo.integrated, species = "astyanax", cluster = cave.names, sum.div = 10, plot.cut = 20, remove.non = T, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600)

plotSankey(object = hypo.integrated, species = "astyanax", cluster = surface.names, sum.div = 10, plot.cut = 20, remove.non = T, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600)


