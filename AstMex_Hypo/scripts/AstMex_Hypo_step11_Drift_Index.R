library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(stringr)
library(viridis)
library(ggpubr)



# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

Idents(hypo.integrated) <- "species.2"

hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

Idents(hypo.integrated) <- "species"

hypo.integrated.surface <- subset(hypo.integrated, idents = "astyanax_surface")
hypo.integrated.cave <- subset(hypo.integrated, idents = "astyanax_cave")

## Load marker gene lists

gene.lists <- readRDS(file = "drift_gene_lists.rds")

names(gene.lists) <- c("conserved.markers", "zebrafish.markers", "astyanax.markers", "conserved.markers.sub", "zebrafish.markers.sub", "astyanax.markers.sub", "conserved.markers.ast", "surface.markers", "cave.markers", "conserved.markers.ast.sub", "surface.markers.sub", "cave.markers.sub")

# ## Need to remove GABA_5 markers, and rename GABA_6 to GABA_5
# ## First remove GABA_5

# for (i in 1:length(gene.lists)) {
	# gene.lists[[i]] <- gene.lists[[i]][!grepl("GABA_5", names(gene.lists[[i]]))]
# }

# for (i in 1:length(gene.lists)) {
	# names(gene.lists[[i]]) <- str_replace(names(gene.lists[[i]]), "GABA_6", "GABA_5")
# }


## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4)] <- lapply(gene.lists[c(1,4)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$zebrafish_avg_logFC > 0 & x[[y]]$astyanax_avg_logFC > 0,]))
gene.lists.pos[c(2,3,5,6,8,9,11,12,14,15,17,18)] <- lapply(gene.lists[c(2,3,5,6,8,9,11,12,14,15,17,18)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))
gene.lists.pos[c(7,10,13,16)] <- lapply(gene.lists[c(7,10,13,16)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}


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


## Plot Subtypes
drift.index.Subtypes <- calcDriftIndex(conserved = gene.lists.pos[[13]], species.1 = gene.lists.pos[[14]], species.2 = gene.lists.pos[[15]])
DI <- tibble(Subtype = c(names(unlist(drift.index.Subtypes))), values = c(unlist(drift.index.Subtypes)))
# ggplot(DI, aes(y = Subtype, x = values, label = Subtype, color = values)) + geom_point(size = 4) + guides(color = FALSE) + ylab(element_blank()) + xlab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

## Plot SubclusterType
drift.index.SubclusterType <- calcDriftIndex(conserved = gene.lists.pos[[16]], species.1 = gene.lists.pos[[17]], species.2 = gene.lists.pos[[18]])
DI.sub <- tibble(SubclusterType = c(names(unlist(drift.index.SubclusterType))), values = c(unlist(drift.index.SubclusterType)))
# ggplot(DI.ast, aes(y = Subtype, x = values, label = Subtype, color = values)) + geom_point(size = 4) + guides(color = FALSE) + ylab(element_blank()) + xlab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + coord_flip()

saveRDS(DI, file = "/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/Drift-index_Subtypes_Amex.rds")
saveRDS(DI.sub, file = "/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/Drift-index_SubclusterTypes_Amex.rds")
 
# OK, make a df with the Subtype as a column
# Plot subtype values ontop!

index <- unique(tibble(Subtype = hypo.integrated.ast@meta.data$Subtype, SubclusterType = hypo.integrated.ast@meta.data$SubclusterType))
DI.sub$Subtype <- index$Subtype[match(DI.sub$SubclusterType, index$SubclusterType)]
DI.sub$Subtype <- factor(DI.sub$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated/Ventricular", "Ependymal", "Progenitors_1", "Progenitors_2", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes_1", "Oligodendrocytes_2", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "GABA_6", "GABA_7", "GABA_Prdx1_1", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Galanin", "Otpa/b_1", "Otpa/b_2", "Otpa/b_3", "Lymphatic", "Tcells", "Bcells", "Mast_cells", "Thrombocytes", "Neutrophils", "Macrophages", "Microglia"))


DI.plot <- ggplot(DI.sub, aes(x = Subtype, y = values, colour = Subtype, fill = Subtype, group = Subtype)) + geom_boxplot(width = 0.5) + geom_point(size = 1, color = "black") + guides(color = FALSE, fill = FALSE) + xlab(element_blank()) + ylab("Drift Index") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DI.plot <- DI.plot + geom_point(data = DI, aes(x = Subtype, y = values), size = 3, shape = 18, color = "red", position = position_nudge(x=0.4)) + theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 8), axis.text.x = element_text(size = 8), axis.title = element_text(size = 10))

png("Figures/Drift_Index_species_subtype.png", width = 8, height = 5, units = "in", res = 250)
DI.plot
dev.off()



# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())][4:9]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files()[grep("GO", list.files())][4:9]
names(go_lists)
go_lists[[7]] <- Reduce(union, list(go_lists[[3]], go_lists[[4]], go_lists[[5]]))


## Plot Subclusters

test <- lapply(go_lists, function(x) c(unlist(calcDriftIndex(conserved = gene.lists.pos[[13]], species.1 = gene.lists.pos[[14]], species.2 = gene.lists.pos[[15]], subset = x, invert = F))))
names(test) <- names(go_lists)

test <- Reduce(cbind, test)
colnames(test) <- names(go_lists)

DI2 <- tibble(cell_type = c(names(unlist(drift.index.SubclusterType))), all = c(unlist(calcDriftIndex(conserved = gene.lists.pos[[16]], species.1 = gene.lists.pos[[17]], species.2 = gene.lists.pos[[18]], invert = F))), TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos[[16]], species.1 = gene.lists.pos[[17]], species.2 = gene.lists.pos[[18]], subset = go_lists[[2]], invert = F))), NP_NTS = c(unlist(calcDriftIndex(conserved = gene.lists.pos[[16]], species.1 = gene.lists.pos[[17]], species.2 = gene.lists.pos[[18]], subset = go_lists[[7]], invert = F))))

DI2 <- DI2[!is.na(DI2$NP_NTS),]

DI3 <- melt(DI2)

ggplot(DI3, aes(x = variable, y = value, group = variable, color = variable)) + geom_jitter(size = 1) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_manual(values = c("blue", "red", "green")) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means()




