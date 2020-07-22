library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

# Load and generate subsets

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

DefaultAssay(hypo.integrated) <- "RNA"

Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

Idents(hypo.integrated) <- "species"
hypo.integrated.surface <- subset(hypo.integrated, idents = "astyanax_surface")
hypo.integrated.cave <- subset(hypo.integrated, idents = "astyanax_cave")

############## Subtype Markers ##############

# Find Conserved Markers

Idents(hypo.integrated.zeb) <- "Subtype"
Idents(hypo.integrated.ast) <- "Subtype"
Idents(hypo.integrated.surface) <- "Subtype"
Idents(hypo.integrated.cave) <- "Subtype"


conserved.markers.ast <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers.ast) <- levels(Idents(hypo.integrated.ast))

# Find Markers for other subsets

surface.markers <- lapply(levels(Idents(hypo.integrated.surface)), function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = F, max.cells.per.ident = 1000))
names(surface.markers) <- levels(Idents(hypo.integrated.surface))

cave.markers <- lapply(levels(Idents(hypo.integrated.cave)), function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = F, max.cells.per.ident = 1000))
names(cave.markers) <- levels(Idents(hypo.integrated.cave))


############## SubclusterType Markers ##############

## Find Conserved Sub Markers
## First identify those subclusters that are species specific

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

Idents(hypo.integrated) <- "SubclusterType"
Idents(hypo.integrated.ast) <- "SubclusterType"
Idents(hypo.integrated.surface) <- "SubclusterType"
Idents(hypo.integrated.cave) <- "SubclusterType"


subclusters <- levels(Idents(hypo.integrated.ast))
index <- vector()
for(i in 1:length(levels(Idents(hypo.integrated.ast)))) {
	scells <- nrow(hypo.integrated.surface@meta.data[hypo.integrated.surface@meta.data$SubclusterType == subclusters[[i]],])
	ccells <- nrow(hypo.integrated.cave@meta.data[hypo.integrated.cave@meta.data$SubclusterType == subclusters[[i]],])
	if (scells < 3 | ccells < 3){
		index[i] <- FALSE
	} else { index[i] <- TRUE}
}

conserved.markers.sub <- lapply(levels(Idents(hypo.integrated.ast))[index], function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = F, max.cells.per.ident = 500))
names(conserved.markers.sub) <- levels(Idents(hypo.integrated.ast))[index]

# Find Sub Markers for each species morph

surface.markers.sub <- lapply(levels(Idents(hypo.integrated.ast))[index], function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(surface.markers.sub) <- levels(Idents(hypo.integrated.ast))[index]

cave.markers.sub <- lapply(levels(Idents(hypo.integrated.ast))[index], function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(cave.markers.sub) <- levels(Idents(hypo.integrated.ast))[index]

# Load gene.lists and append new marker gene lists

gene.lists <- readRDS("drift_gene_lists.rds")


conserved.markers.ast.subs <- conserved.markers.ast
surface.markers.ast.subs <- surface.markers
cave.markers.ast.sub <- cave.markers
conserved.markers.sub.ast.sub <- conserved.markers.sub
surface.markers.sub.ast.sub <- surface.markers.sub
cave.markers.sub.ast.sub <- cave.markers.sub

new.lists <- list(conserved.markers.ast.subs, surface.markers.ast.subs, cave.markers.ast.sub, conserved.markers.sub.ast.sub, surface.markers.sub.ast.sub, cave.markers.sub.ast.sub)
names(new.lists) <- c("conserved.markers.ast.subs", "surface.markers.ast.subs", "cave.markers.ast.sub", "conserved.markers.sub.ast.sub", "surface.markers.sub.ast.sub", "cave.markers.sub.ast.sub")

gene.lists <- c(gene.lists, new.lists)

saveRDS(gene.lists, file = "drift_gene_lists.rds")
print("done cave markers")






