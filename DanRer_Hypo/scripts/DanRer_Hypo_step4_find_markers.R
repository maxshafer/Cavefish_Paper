library(Seurat)
library(Matrix)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

## Load objects

hypo <- readRDS("DanRer_65k.rds")

############## Subtype Markers ##############

# Find Conserved Markers

Idents(hypo) <- "Subtype"

subtype.markers <- lapply(levels(Idents(hypo)), function(x) FindMarkers(hypo, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(subtype.markers) <- levels(Idents(hypo))

gene.lists <- list()

gene.lists[[1]] <- subtype.markers

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("Done Subtype marker gene lists")

############## SubclusterType Markers ##############

## Find SubcluserType markers

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

gene.lists <- readRDS("marker_gene_lists.rds")

Idents(hypo) <- "SubclusterType"

subcluster.markers <- lapply(Idents(hypo), function(x) FindMarkers(hypo, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(subcluster.markers) <- index

gene.lists[[2]] <- subcluster.markers

saveRDS(gene.lists, file = "marker_gene_lists.rds")
print("done SubclusterType markers")



