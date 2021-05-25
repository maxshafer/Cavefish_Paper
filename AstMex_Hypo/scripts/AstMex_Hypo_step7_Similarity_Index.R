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

## Calculate the similarity index!
## Takes 3 lists of marker genes, conserved, speices.1 and species.2
## Takes the intersection of all 3 Cluster names, and subset/order lists by that
## Then run Similarity Index caculation for each Cluster
## Maybe add a part for using a cutoff (p-val, logfc, percent)
## Add a part for putting number of genes

### Should I add a part where I make other columns in SI for number of marker genes? Should be easy (relatively)

calcSimilarityIndex <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	SI <- list()
	if (is.null(subset)) {
		SI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (nrow(conserved[[x]]) / nrow(species.1[[x]]))) * (1 - (nrow(conserved[[x]]) / nrow(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		SI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		SI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
	}
	}
	names(SI) <- names
	return(SI)
}


## SI Clusters
similarity.index.Clusters <- calcSimilarityIndex(conserved = gene.lists.pos[[1]], species.1 = gene.lists.pos[[2]], species.2 = gene.lists.pos[[3]])
SI <- tibble(Cluster = c(names(unlist(similarity.index.Clusters))), values = c(unlist(similarity.index.Clusters)))

## SI Subcluster
similarity.index.Subcluster <- calcSimilarityIndex(conserved = gene.lists.pos[[4]], species.1 = gene.lists.pos[[5]], species.2 = gene.lists.pos[[6]])
SI.sub <- tibble(Subcluster = c(names(unlist(similarity.index.Subcluster))), values = c(unlist(similarity.index.Subcluster)))

saveRDS(SI, file = "Ast_SI.rds")
saveRDS(SI.sub, file = "Ast_SI.sub.rds")
