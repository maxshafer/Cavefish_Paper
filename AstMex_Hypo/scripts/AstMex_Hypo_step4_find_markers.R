library(Seurat)
library(Matrix)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k_vR.rds")

# Subset to cave and surface versions, and find conserved, and species-morph specific marker genes

DefaultAssay(hypo) <- "RNA"

Idents(hypo) <- "species"
hypo.surface <- subset(hypo, idents = "astyanax_surface")
hypo.cave <- subset(hypo, idents = "astyanax_cave")

############## Cluster Markers ##############

# Find Conserved Markers

Idents(hypo) <- "Cluster"
Idents(hypo.surface) <- "Cluster"
Idents(hypo.cave) <- "Cluster"

print("Done loading, starting FindConservedMarkers")

conserved.markers <- lapply(levels(Idents(hypo)), function(x) FindConservedMarkers(hypo, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers) <- levels(Idents(hypo))

gene.lists <- list()

gene.lists[[1]] <- conserved.markers

# Find Markers for species morphs

surface.markers <- lapply(levels(Idents(hypo.surface)), function(x) FindMarkers(hypo.surface, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(surface.markers) <- levels(Idents(hypo.surface))

gene.lists[[2]] <- surface.markers

cave.markers <- lapply(levels(Idents(hypo.cave)), function(x) FindMarkers(hypo.cave, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(cave.markers) <- levels(Idents(hypo.cave))

gene.lists[[3]] <- cave.markers

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("Done Cluster marker gene lists")

############## Subcluster Markers ##############

## Find Conserved Sub Markers
## First identify those subclusters that are specific (> 90%), I don't care about DE genes for these cells types

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

gene.lists <- readRDS("marker_gene_lists.rds")

Idents(hypo) <- "Subcluster"
Idents(hypo.surface) <- "Subcluster"
Idents(hypo.cave) <- "Subcluster"


# Find morph specific clusters
prop.table <- table(hypo@meta.data$Subcluster, hypo@meta.data$species)
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo))
index <- subclusters[!(subclusters %in% c(surface.names, cave.names))]

## Find conserved markers for shared cell types, and markers for specific cell types

print("starting conserved markers")

conserved.markers.sub <- lapply(index, function(x) FindConservedMarkers(hypo, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 500))
names(conserved.markers.sub) <- index

conserved.markers.sub.specific <- lapply(c(surface.names, cave.names), function(x) FindMarkers(hypo, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(conserved.markers.sub.specific) <- c(surface.names, cave.names)

conserved.markers.sub <- c(conserved.markers.sub, conserved.markers.sub.specific)
conserved.markers.sub <- conserved.markers.sub[levels(Idents(hypo))]

gene.lists[[4]] <- conserved.markers.sub

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done conserved Subcluster markers")

# Find markers for surface and cave versions of each subcluster

surface.markers.sub <- lapply(index, function(x) FindMarkers(hypo.surface, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(surface.markers.sub) <- index

surface.markers.sub.specific <- lapply(surface.names, function(x) FindMarkers(hypo.surface, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(surface.markers.sub.specific) <- surface.names

surface.markers.sub <- c(surface.markers.sub, surface.markers.sub.specific)
surface.markers.sub <- surface.markers.sub[levels(Idents(hypo))]

gene.lists[[5]] <- surface.markers.sub

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done surface Subcluster markers")


cave.markers.sub <- lapply(index, function(x) FindMarkers(hypo.cave, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(cave.markers.sub) <- index

cave.markers.sub.specific <- lapply(cave.names, function(x) FindMarkers(hypo.cave, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(cave.markers.sub.specific) <- cave.names

cave.markers.sub <- c(cave.markers.sub, cave.markers.sub.specific)
cave.markers.sub <- cave.markers.sub[levels(Idents(hypo))]

gene.lists[[6]] <- cave.markers.sub

names(gene.lists) <- c("cluster.conserved", "cluster.surface", "cluster.cave", "subcluster.conserved", "subcluster.surface", "subcluster.cave")

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done cave Subcluster markers")

# Find markers for surface and cave versions of each subcluster

Idents(hypo) <- "Cluster"
de.markers <- lapply(levels(Idents(hypo)), function(x) {
  subset <- subset(hypo, ident = x)
  Idents(subset) <- "species"
  markers <- FindMarkers(subset, ident.1 = "astyanax_surface", ident.2 = "astyanax_cave", verbose = T, max.cells.per.ident = 500)
  return(markers)
})
names(de.markers) <- levels(Idents(hypo))

gene.lists[[7]] <- de.markers

Idents(hypo) <- "Subcluster"
de.markers.sub <- lapply(index, function(x) {
  subset <- subset(hypo, ident = x)
  Idents(subset) <- "species"
  markers <- FindMarkers(subset, ident.1 = "astyanax_surface", ident.2 = "astyanax_cave", verbose = T, max.cells.per.ident = 500)
  return(markers)
})
names(de.markers.sub) <- index

gene.lists[[8]] <- de.markers.sub

names(gene.lists)[7:8] <- c("cluster.de", "subcluster.de")
saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done DE Subcluster markers")
