library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggpubr)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k_vR.rds")

Idents(hypo) <- "species"

hypo.surface <- subset(hypo, idents = "astyanax_surface")
hypo.cave <- subset(hypo, idents = "astyanax_cave")

## Load marker gene lists

gene.lists <- readRDS(file = "marker_gene_lists.rds")

## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4)] <- lapply(gene.lists[c(1,4)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))
gene.lists.pos[c(2,3,5,6)] <- lapply(gene.lists[c(2,3,5,6)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)[c(1:6)]

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

saveRDS(gene.lists.pos, file = "marker_gene_lists_pos.rds")

## Calculate the drift index!
## Takes 3 lists of marker genes, conserved, speices.1 and species.2
## Takes the intersection of all 3 Cluster names, and subset/order lists by that
## Then run Drift Index caculation for each Cluster
## Maybe add a part for using a cutoff (p-val, logfc, percent)
## Add a part for putting number of genes

### Should I add a part where I make other columns in DI for number of marker genes? Should be easy (relatively)

calcDriftIndex <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	DI <- list()
	if (is.null(subset)) {
		DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (nrow(conserved[[x]]) / nrow(species.1[[x]]))) * (1 - (nrow(conserved[[x]]) / nrow(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
	}
	}
	names(DI) <- names
	return(DI)
}


## DI Clusters
drift.index.Clusters <- calcDriftIndex(conserved = gene.lists.pos[[1]], species.1 = gene.lists.pos[[2]], species.2 = gene.lists.pos[[3]])
DI <- tibble(Cluster = c(names(unlist(drift.index.Clusters))), values = c(unlist(drift.index.Clusters)))

## DI Subcluster
drift.index.Subcluster <- calcDriftIndex(conserved = gene.lists.pos[[4]], species.1 = gene.lists.pos[[5]], species.2 = gene.lists.pos[[6]])
DI.sub <- tibble(Subcluster = c(names(unlist(drift.index.Subcluster))), values = c(unlist(drift.index.Subcluster)))

saveRDS(DI, file = "Ast_DI.rds")
saveRDS(DI.sub, file = "Ast_DI.sub.rds")
