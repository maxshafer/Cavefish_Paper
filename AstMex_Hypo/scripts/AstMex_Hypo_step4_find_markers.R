library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k.rds")

# Subset to cave and surface versions, and find conserved, and species-morph specific marker genes

DefaultAssay(hypo) <- "RNA"

Idents(hypo) <- "species"
hypo.surface <- subset(hypo, idents = "astyanax_surface")
hypo.cave <- subset(hypo, idents = "astyanax_cave")

# ############## Subtype Markers ##############
# 
# # Find Conserved Markers
# 
# Idents(hypo) <- "Subtype"
# Idents(hypo.surface) <- "Subtype"
# Idents(hypo.cave) <- "Subtype"
# 
# conserved.markers <- lapply(levels(Idents(hypo)), function(x) FindConservedMarkers(hypo, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 1000))
# names(conserved.markers) <- levels(Idents(hypo))
# 
# gene.lists <- list()
# 
# gene.lists[[1]] <- conserved.markers
# 
# # Find Markers for species morphs
# 
# surface.markers <- lapply(levels(Idents(hypo.surface)), function(x) FindMarkers(hypo.surface, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
# names(surface.markers) <- levels(Idents(hypo.surface))
# 
# gene.lists[[2]] <- surface.markers
# 
# cave.markers <- lapply(levels(Idents(hypo.cave)), function(x) FindMarkers(hypo.cave, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
# names(cave.markers) <- levels(Idents(hypo.cave))
# 
# gene.lists[[3]] <- cave.markers
# 
# names(gene.lists) <- c("conserved.markers", "surface.markers", "cave.markers")
# 
# saveRDS(gene.lists, file = "marker_gene_lists.rds")

############## SubclusterType Markers ##############

## Find Conserved Sub Markers
## First identify those subclusters that are specific (> 90%), I don't care about DE genes for these cells types

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

Idents(hypo) <- "SubclusterType"
Idents(hypo.surface) <- "SubclusterType"
Idents(hypo.cave) <- "SubclusterType"


# Find morph specific clusters
prop.table <- table(hypo@meta.data$SubclusterType, hypo@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$surface_specific <- ifelse(prop.table$astyanax_surface > .9, "yes", "no")
prop.table$cave_specific <- ifelse(prop.table$astyanax_cave > .9, "yes", "no")

surface.names <- row.names(prop.table[prop.table$surface_specific == "yes",])
cave.names <- row.names(prop.table[prop.table$cave_specific == "yes",])

# Make index
subclusters <- levels(Idents(hypo))
index <- subclusters[!(subclusters %in% c(surface.names, cave.names))]

## Find conserved markers for shared cell types, and markers for specific cell types

print("starting conserved markers")

conserved.markers.sub <- lapply(index, function(x) FindConservedMarkers(hypo, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 500))
names(conserved.markers.sub) <- index

conserved.markers.sub.specific <- lapply(c(surface.names, cave.names), function(x) FindMarkers(hypo, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(conserved.markers.sub.specific) <- c(surface.names, cave.names)

conserved.markers.sub <- c(conserved.markers.sub, conserved.markers.sub.specific)
conserved.markers.sub <- conserved.markers.sub[levels(Idents(hypo))]

gene.lists[[4]] <- conserved.markers.sub

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done conserved markers")

# Find markers for surface and cave versions of each subcluster

surface.markers.sub <- lapply(levels(Idents(hypo))[index], function(x) FindMarkers(hypo.surface, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(surface.markers.sub) <- levels(Idents(hypo))[index]

surface.markers.sub.specific <- lapply(surface.names, function(x) FindMarkers(hypo.surface, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(surface.markers.sub.specific) <- surface.names

surface.markers.sub <- c(surface.markers.sub, surface.markers.sub.specific)
surface.markers.sub <- surface.markers.sub[levels(Idents(hypo))]

gene.lists[[5]] <- surface.markers.sub

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done surface markers")


cave.markers.sub <- lapply(levels(Idents(hypo))[index], function(x) FindMarkers(hypo.cave, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(cave.markers.sub) <- levels(Idents(hypo))[index]

cave.markers.sub.specific <- lapply(cave.names, function(x) FindMarkers(hypo.cave, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(cave.markers.sub.specific) <- cave.names

cave.markers.sub <- c(cave.markers.sub, cave.markers.sub.specific)
cave.markers.sub <- cave.markers.sub[levels(Idents(hypo))]

gene.lists[[6]] <- cave.markers.sub

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done cave markers")



