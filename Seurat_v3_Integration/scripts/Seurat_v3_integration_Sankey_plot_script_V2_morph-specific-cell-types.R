library(networkD3)
library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v2.rds")

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?
Idents(hypo.integrated) <- "species.2"
hypo.ast <- subset(hypo.integrated, idents = "astyanax")

prop.table <- table(hypo.ast@meta.data$SubclusterType, hypo.ast@meta.data$species)
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Function to plot Sankey's for each
## object : the seurat object used to extract cell labels, defaults to "hypo.integrated"
## species : defaults to "asytanax"
## cluster : accepts single character value, or a vector or characters indicating the species-morph specific clusters (from Subtype, or SubclusterType), not integrated idents
## sum.div : cutoff for which integrated identies to include - a value of 10 will include only those integrated identities where at least 10% of the "cluster" cells map. 20 will be 5%, and 5 will be 20%
## plot.cut : minimum number of cells per connection for plotting
## remove.non : Binary. Set to TRUE to remove astyanax clusters that are not in "cluster" above. If FALSE, astyanax clusters that also map to the integrated idents that "cluster" maps to will be included
## check.species.2 : accepts "yes" or "no", if "yes", will check if the number of cells connecting to each zebrafish identity is greater than 10% of the total number of zebrafish cells in that identity

plotSankey <- function(object = hypo.integrated, species = "astyanax", cluster = "Glut_1_4", check.species.2 = "no", sum.div = 10, plot.cut = 10, remove.non = F, ...) {
  
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
  object2 <- subset(object, idents = subset.names)
  
  # Extract the count tables
  if (res == "SubclusterType") {
    if (!(species %in% object2@meta.data$species.2)) {
      table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, "integrated_SubclusterType"], object2@meta.data[object2@meta.data$species.2 != species,"SubclusterType"]))
      table1 <- table2
      table1$value <- 0
      print(paste(cluster, "does not exist in", species, sep = " "))
    } else {
      if (!(species.2 %in% object2@meta.data$species.2)) {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, "SubclusterType"], object2@meta.data[object2@meta.data$species.2 == species,"integrated_SubclusterType"]))
        table2 <- table1
        table2$value <- 0
        print(paste(cluster, "does not exist in", species.2, sep = " "))
      } else {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, "SubclusterType"], object2@meta.data[object2@meta.data$species.2 == species,"integrated_SubclusterType"]))
        table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, "integrated_SubclusterType"], object2@meta.data[object2@meta.data$species.2 != species,"SubclusterType"]))
      }
    }
  }
  if (res == "Subtype") {
    if (!(species %in% object2@meta.data$species.2)) {
      table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, "integrated_Subtype"], object2@meta.data[object2@meta.data$species.2 != species,"Subtype"]))
      table1 <- table2
      table1$value <- 0
      print(paste(cluster, "does not exist in", species, sep = " "))
    } else {
      if (!(species.2 %in% object2@meta.data$species.2)) {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, "Subtype"], object2@meta.data[object2@meta.data$species.2 == species,"integrated_Subtype"]))
        table2 <- table1
        table2$value <- 0
        print(paste(cluster, "does not exist in", species.2, sep = " "))
      } else {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, "Subtype"], object2@meta.data[object2@meta.data$species.2 == species,"integrated_Subtype"]))
        table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, "integrated_Subtype"], object2@meta.data[object2@meta.data$species.2 != species,"Subtype"]))
      }
    }
  }
  
  colnames(table1) <- c("Var.1", "Var.2", "value")
  colnames(table2) <- c("Var.1", "Var.2", "value")
  
  # Check whether the zebrafish cluster connections represent at least 10% of that cluster
  if (check.species.2 == "yes") {
    species.2.table <- table(object@meta.data[object@meta.data$species.2 != species, "SubclusterType"])
    index <- table2$value > species.2.table[table2$Var.2]/10
    table2 <- table2[index,]
  }
  
  # Append species idents + combine lists + remove lowbounds
  table1$Var1 <- paste(as.character(species), table1$Var.1, sep = "_")
  table1$Var2 <- paste("int", table1$Var.2, sep = "_")
  table2$Var1 <- paste("int", table2$Var.1, sep = "_")
  table2$Var2 <- paste(as.character(species.2), table2$Var.2, sep = "_")
  
  table <- rbind(table1, table2)
  table <- table[table$value > plot.cut,]
  
  if (remove.non == TRUE) {
    cluster2 <- paste(species, cluster, sep = "_")
    table <- table[c(seq(1:length(table$Var1))[table$Var1 %in% cluster2], grep("int", table$Var1)),]
  }
  
  # Make nodes list
  nodes <- union(as.character(table$Var2), as.character(table$Var1))
  nodes <- as.data.frame(nodes)
  
  # Fix node levels
  table$Var.1 <- as.numeric(row.names(nodes)[match(table$Var1, nodes$nodes)])-1
  table$Var.2 <- as.numeric(row.names(nodes)[match(table$Var2, nodes$nodes)])-1

  
  return(sankeyNetwork(Links = table, Nodes = nodes, Source = "Var.1", Target = "Var.2", Value = "value", NodeID = "nodes", units = "cells", iterations = 10000, sinksRight = F, ...))
}

## Make some figures
# sum.div = 10, plot.cut = 20 were defaults

# For Figure 6e
plotSankey(object = hypo.integrated, species = "astyanax", cluster = "Glut_1_4", check.species.2 = "yes", sum.div = 10, plot.cut = 5, remove.non = F, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 175, width = 400)

# For Figure S8c and S8d
plotSankey(object = hypo.integrated, species = "astyanax", cluster = cave.names, check.species.2 = "yes", sum.div = 5, plot.cut = 25, remove.non = T, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600)

plotSankey(object = hypo.integrated, species = "astyanax", cluster = surface.names, check.species.2 = "yes", sum.div = 10, plot.cut = 25, remove.non = T, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600)


