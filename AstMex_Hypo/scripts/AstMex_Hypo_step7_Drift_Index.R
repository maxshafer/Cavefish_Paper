library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)
library(ggpubr)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k.rds")

Idents(hypo) <- "species"

hypo.surface <- subset(hypo, idents = "astyanax_surface")
hypo.cave <- subset(hypo, idents = "astyanax_cave")

## Load marker gene lists

gene.lists <- readRDS(file = "marker_gene_lists.rds")

## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4)] <- lapply(gene.lists[c(1,4)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))
gene.lists.pos[c(2,3,5,6)] <- lapply(gene.lists[c(2,3,5,6)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

saveRDS(gene.lists.pos, file = "marker_gene_lists_pos.rds")

## Calculate the drift index!
## Takes 3 lists of marker genes, conserved, speices.1 and species.2
## Takes the intersection of all 3 Subtype names, and subset/order lists by that
## Then run Drift Index caculation for each Subtype
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
		DI <- lapply(seq_along(conserved), function(x) sqrt( abs( (1 - (nrow(conserved[[x]]) / nrow(species.1[[x]]))) * (1 - (nrow(conserved[[x]]) / nrow(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		DI <- lapply(seq_along(conserved), function(x) sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		DI <- lapply(seq_along(conserved), function(x) sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
	}
	}
	names(DI) <- names
	return(DI)
}


## DI Subtypes
drift.index.Subtypes <- calcDriftIndex(conserved = gene.lists.pos[[1]], species.1 = gene.lists.pos[[2]], species.2 = gene.lists.pos[[3]])
DI <- tibble(Subtype = c(names(unlist(drift.index.Subtypes))), values = c(unlist(drift.index.Subtypes)))

## DI SubclusterType
drift.index.SubclusterType <- calcDriftIndex(conserved = gene.lists.pos[[4]], species.1 = gene.lists.pos[[5]], species.2 = gene.lists.pos[[6]])
DI.sub <- tibble(SubclusterType = c(names(unlist(drift.index.SubclusterType))), values = c(unlist(drift.index.SubclusterType)))

saveRDS(DI, file = "Ast_DI.rds")
saveRDS(DI.sub, file = "Ast_DI.sub.rds")
 
# OK, make a df with the Subtype as a column
# Plot subtype values ontop!

index <- unique(tibble(Subtype = hypo@meta.data$Subtype, SubclusterType = hypo@meta.data$SubclusterType))
DI.sub$Subtype <- index$Subtype[match(DI.sub$SubclusterType, index$SubclusterType)]
DI.sub$Subtype <- factor(DI.sub$Subtype, levels = levels(hypo@meta.data$Subtype))


DI.plot <- ggplot(DI.sub, aes(x = Subtype, y = values, colour = Subtype, fill = Subtype, group = Subtype)) + geom_boxplot(width = 0.5) + geom_point(size = 1, color = "black") + guides(color = FALSE, fill = FALSE) + xlab(element_blank()) + ylab("Drift Index") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DI.plot <- DI.plot + geom_point(data = DI, aes(x = Subtype, y = values), size = 3, shape = 18, color = "red", position = position_nudge(x=0.4)) + theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 8), axis.text.x = element_text(size = 8), axis.title = element_text(size = 10))

png("Figures/Drift_Index_species_subtype.png", width = 8, height = 5, units = "in", res = 250)
DI.plot
dev.off()



# Get Go lists and use only those that are marker genes

go_lists <- list.files("../GO_lists/")[grep("GO", list.files("../GO_lists/"))]
go_lists <- lapply(go_lists, function(x) read.csv(paste("../GO_lists/", x, sep = ""), head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files("../GO_lists/")[grep("GO", list.files("../GO_lists/"))]
names(go_lists)
go_lists[[7]] <- Reduce(union, list(go_lists[[3]], go_lists[[4]], go_lists[[5]]))


## Plot Subclusters

test <- lapply(go_lists, function(x) c(unlist(calcDriftIndex(conserved = gene.lists.pos[[1]], species.1 = gene.lists.pos[[2]], species.2 = gene.lists.pos[[3]], subset = x, invert = F))))
names(test) <- names(go_lists)

test <- Reduce(cbind, test)
colnames(test) <- names(go_lists)

DI2 <- tibble(cell_type = c(names(unlist(drift.index.SubclusterType))), all = c(unlist(calcDriftIndex(conserved = gene.lists.pos[[4]], species.1 = gene.lists.pos[[5]], species.2 = gene.lists.pos[[6]], invert = F))), TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos[[4]], species.1 = gene.lists.pos[[5]], species.2 = gene.lists.pos[[6]], subset = go_lists[[2]], invert = F))), NP_NTS = c(unlist(calcDriftIndex(conserved = gene.lists.pos[[4]], species.1 = gene.lists.pos[[5]], species.2 = gene.lists.pos[[6]], subset = go_lists[[7]], invert = F))))
DI2 <- DI2[!is.na(DI2$NP_NTS),]
DI3 <- melt(DI2)

ggplot(DI3, aes(x = variable, y = value, group = variable, color = variable)) + geom_jitter(size = 1) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_manual(values = c("blue", "red", "green")) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means()
