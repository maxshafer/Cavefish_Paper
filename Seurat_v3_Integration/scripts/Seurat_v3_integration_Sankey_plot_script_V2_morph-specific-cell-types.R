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
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

# Find species specific clusters
prop.table <- table(hypo.integrated@meta.data$integrated_Subcluster, hypo.integrated@meta.data$species.2)
zeb.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
ast.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

prop.table <- table(hypo.integrated.ast@meta.data$Subcluster, hypo.integrated.ast@meta.data$species)
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Function to plot Sankey's for each
## object : the seurat object used to extract cell labels, defaults to "hypo.integrated"
## species : defaults to "asytanax"
## integrated : TRUE or FALSE, determines if it should look in the integrated identities (or not) for the provided 'cluster'(s)
## cluster : accepts single character value, or a vector or characters indicating the species-morph specific clusters (from Cluster, or Subcluster), not integrated idents
## sum.div : cutoff for which integrated identies to include - a value of 10 will include only those integrated identities where at least 10% of the "cluster" cells map. 20 will be 5%, and 5 will be 20%
## plot.cut : minimum number of cells per connection for plotting
## remove.non : Binary. Set to TRUE to remove astyanax clusters that are not in "cluster" above. If FALSE, astyanax clusters that also map to the integrated idents that "cluster" maps to will be included
## check.species.2 : accepts "yes" or "no", if "yes", will check if the number of cells connecting to each zebrafish identity is greater than 10% of the total number of zebrafish cells in that identity

plotSankey <- function(object = hypo.integrated, species = "astyanax", integrated = FALSE, cluster = "Neuronal_03_4", check.species.2 = "no", sum.div = 10, plot.cut = 10, remove.non = F, ...) {
  
  if (species == "astyanax") {
    species.2 <- "zebrafish"
  } else {
    species.2 <- "astyanax"
  }
  
  # Set source and targets (if using integrated idents as input)
  if (integrated == FALSE) {
    # using independent idents as input
    source.cluster <- "Cluster"
    target.cluster <- "integrated_Cluster"
    source.subcluster <- "Subcluster"
    target.subcluster <- "integrated_Subcluster"
  } 
  # else {
  #   # using integrated idents as input
  #   source.cluster <- "integrated_Cluster"
  #   target.cluster <- "Cluster"
  #   source.subcluster <- "integrated_Subcluster"
  #   target.subcluster <- "Subcluster"
  # }
  
  # Identify the integrated IDs that the target cluster maps to
  if (integrated == FALSE) {
    if (cluster %in% unique(object@meta.data[,source.subcluster])) {
      res <- source.subcluster
      Idents(object) <- target.subcluster
      if (length(cluster) == 1) {
        table <- table(object@meta.data[object@meta.data$species.2 == species & object@meta.data[,source.subcluster] %in% cluster, target.subcluster])
        subset.names <- names(table[table > sum(table)/sum.div])
      }
      if (length(cluster) > 1) {
        subset.names <- list()
        for (i in 1:length(cluster)) {
          table <- table(object@meta.data[object@meta.data$species.2 == species & object@meta.data[,source.subcluster] %in% cluster[[i]], target.subcluster])
          subset.names[[i]] <- names(table[table > sum(table)/sum.div])
        }
        subset.names <- unlist(subset.names)
      }
    }
  }
   
  
  if (cluster %in% unique(object@meta.data[,source.cluster])) {
    res <- source.cluster
    Idents(object) <- target.cluster
    table <- table(object@meta.data[object@meta.data$species.2 == species & object@meta.data[,source.cluster] %in% cluster, source.cluster])
  }
  
  # Subset object to mapped integrated IDs
  #subset.names <- names(table[table > sum(table)/20]) # The problem: if using multiple idents, this should be normalized for the size of each
  if (integrated == FALSE) {
    print("Integrated identities used for plotting")
  } else {
    print(paste(species, "identities used for plotting", sep = " "))
    subset.names <- cluster
    Idents(object) <- "integrated_Subcluster"
  }
  print(subset.names)
  object2 <- subset(object, idents = subset.names)
  
  # Extract the count tables
  if (res == "Subcluster" | res == "integrated_Subcluster") {
    if (!(species %in% object2@meta.data$species.2)) {
      table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, target.subcluster], object2@meta.data[object2@meta.data$species.2 != species,source.subcluster]))
      table1 <- table2
      table1$value <- 0
      print(paste(cluster, "does not exist in", species, sep = " "))
    } else {
      if (!(species.2 %in% object2@meta.data$species.2)) {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, source.subcluster], object2@meta.data[object2@meta.data$species.2 == species,target.subcluster]))
        table2 <- table1
        table2$value <- 0
        print(paste(cluster, "does not exist in", species.2, sep = " "))
      } else {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, source.subcluster], object2@meta.data[object2@meta.data$species.2 == species,target.subcluster]))
        table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, target.subcluster], object2@meta.data[object2@meta.data$species.2 != species,source.subcluster]))
      }
    }
  }
  if (res == "Cluster" | res == "integrated_Cluster") {
    if (!(species %in% object2@meta.data$species.2)) {
      table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, target.cluster], object2@meta.data[object2@meta.data$species.2 != species,source.cluster]))
      table1 <- table2
      table1$value <- 0
      print(paste(cluster, "does not exist in", species, sep = " "))
    } else {
      if (!(species.2 %in% object2@meta.data$species.2)) {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, source.cluster], object2@meta.data[object2@meta.data$species.2 == species,target.cluster]))
        table2 <- table1
        table2$value <- 0
        print(paste(cluster, "does not exist in", species.2, sep = " "))
      } else {
        table1 <- melt(table(object2@meta.data[object2@meta.data$species.2 == species, source.cluster], object2@meta.data[object2@meta.data$species.2 == species,target.cluster]))
        table2 <- melt(table(object2@meta.data[object2@meta.data$species.2 != species, target.cluster], object2@meta.data[object2@meta.data$species.2 != species,source.cluster]))
      }
    }
  }
  
  colnames(table1) <- c("Var.1", "Var.2", "value")
  colnames(table2) <- c("Var.1", "Var.2", "value")
  
  # Check whether the zebrafish cluster connections represent at least 10% of that cluster
  if (check.species.2 == "yes") {
    species.2.table <- table(object@meta.data[object@meta.data$species.2 != species, source.subcluster])
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

  
  return(sankeyNetwork(Links = table, Nodes = nodes, Source = "Var.1", Target = "Var.2", Value = "value", NodeID = "nodes", units = "cells", iterations = 10000, ...))
}

## Make some figures
# sum.div = 10, plot.cut = 20 were defaults

# For Figure 6e
plotSankey(object = hypo.integrated, species = "astyanax", cluster = "Neuronal_03_4", check.species.2 = "yes", sum.div = 10, plot.cut = 5, remove.non = F, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 175, width = 400, sinksRight = F)

# For Figure S8c and S8d
plotSankey(object = hypo.integrated, species = "astyanax", cluster = cave.names, check.species.2 = "yes", sum.div = 5, plot.cut = 25, remove.non = T, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600, sinksRight = F)

# For zebrafish specific integrated subclusters (map them back to zebrafish independent subclusters)
plotSankey(object = hypo.integrated, species = "zebrafish", integrated = TRUE, cluster = zeb.names, check.species.2 = "yes", sum.div = 10, plot.cut = 25, remove.non = F, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600, sinksRight = F)


# Make Sankeys + DotPlots for removed subclusters
# First load old data, rename meta data columns so function works
# using no cutoffs for Sankey to demonstrae the complete weirdness

hypo.integrated.v1 <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_v1.rds")
Idents(hypo.integrated.v1) <- "integrated_SubclusterType"

subset <- subset(hypo.integrated.v1, cells = WhichCells(hypo.integrated.v1, idents = c("Glut_2_5", "Glut_2_4", "Glut_3_6", "Glut_5_5", "Glut_3_7", "GABA_1_12")))

colnames(hypo.integrated.v1@meta.data)[c(7,8,21,23)] <- c("Cluster", "Subcluster", "integrated_Cluster", "integrated_Subcluster")
removed.idents <- c("Glut_2_5", "Glut_2_4", "Glut_3_6", "Glut_5_5", "Glut_3_7", "GABA_1_12")

plotSankey(object = hypo.integrated.v1, species = "zebrafish", integrated = T, cluster = removed.idents, check.species.2 = "yes", sum.div = 1, plot.cut = 1, remove.non = F, fontSize = 12, nodeWidth = 40, nodePadding = 10, height = 300, width = 600, sinksRight = T)

subset.markers <- lapply(removed.idents, function(x) FindMarkers(hypo.integrated.v1, ident.1 = x, min.pct = 0.1, max.cells.per.ident = 500))
genes.to.plot <- lapply(subset.markers, function(x) row.names(x)[1:10])

removed.dots <- DotPlot(subset, features = unique(unlist(genes.to.plot)), group.by = "integrated_SubclusterType", scale.max = 200) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title = element_blank())

pdf("Figures/Removed-idents_dot-plot.pdf", width = 10, height = 10)
removed.dots + plot_layout(guides = "collect", height = unit(c(50), c("mm")), width = unit(c(100), c("mm")))
dev.off()




