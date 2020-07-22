library(Seurat)
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

# names(gene.lists) <- c("conserved.markers", "zebrafish.markers", "astyanax.markers", "conserved.markers.sub", "zebrafish.markers.sub", "astyanax.markers.sub", "conserved.markers.ast", "surface.markers", "cave.markers", "conserved.markers.ast.sub", "surface.markers.sub", "cave.markers.sub")

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
gene.lists.pos[c(2,3,5,6,8,9,11,12)] <- lapply(gene.lists[c(2,3,5,6,8,9,11,12)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))
gene.lists.pos[c(7,10)] <- lapply(gene.lists[c(7,10)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists[1:12])

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

## Plot the # of conserved and species specific markers

marker.numbers <- data.frame(subtype = names(unlist(lapply(gene.lists.pos[[1]], function(x) nrow(x)))), conserved = unlist(lapply(gene.lists.pos[[1]], function(x) nrow(x))), zebrafish = unlist(lapply(gene.lists.pos[[2]], function(x) nrow(x))), cavefish = unlist(lapply(gene.lists.pos[[3]], function(x) nrow(x))))

marker.numbers.sub <- data.frame(subcluster = names(unlist(lapply(gene.lists.pos[[4]], function(x) nrow(x)))), subtype = hypo.integrated@meta.data$integrated_Subtype[match(names(unlist(lapply(gene.lists.pos[[4]], function(x) nrow(x)))), hypo.integrated@meta.data$integrated_SubclusterType)], conserved = unlist(lapply(gene.lists.pos[[4]], function(x) nrow(x))), zebrafish = unlist(lapply(gene.lists.pos[[5]], function(x) nrow(x))), cavefish = unlist(lapply(gene.lists.pos[[6]], function(x) nrow(x))))

marker.numbers$subtype <- factor(marker.numbers$subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))
marker.numbers.sub$subtype <- factor(marker.numbers.sub$subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

ggplot(marker.numbers.sub, aes(x = subtype, y = conserved, fill = subtype)) + geom_boxplot(outlier.color = "transparent") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + guides(fill = F) + coord_flip() + geom_jitter(data = marker.numbers.sub, aes(x = subtype, y = conserved), color = "black") + geom_point(data = marker.numbers, aes(x = subtype, y = conserved), color = "red")

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


## Plot the number of conserved and species specific marker genes per cluster

test <- lapply(seq_along(1:length(gene.lists.pos[[4]])), function(x) {
	zeb <- length(setdiff(rownames(gene.lists.pos[[5]][[x]]), rownames(gene.lists.pos[[4]][[x]])))
	ast <- length(setdiff(rownames(gene.lists.pos[[6]][[x]]), rownames(gene.lists.pos[[4]][[x]])))
	con <- length(rownames(gene.lists.pos[[4]][[x]]))
	return(data.frame(con, ast, zeb)) } )


test <- Reduce(rbind, test)
test$cell_type <- names(gene.lists.pos[[4]])

ggplot(melt(test), aes(y = variable, x = cell_type,size = value, color = value)) + geom_point() + scale_color_viridis(option = "B") + coord_flip() + theme(axis.text = element_text(size = 6))

## Plot Subtypes
drift.index.Subtypes <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers, species.1 = gene.lists.pos$zebrafish.markers, species.2 = gene.lists.pos$astyanax.markers)
DI <- tibble(Subtype = c(names(unlist(drift.index.Subtypes))), values = c(unlist(drift.index.Subtypes)))
# ggplot(DI, aes(y = Subtype, x = values, label = Subtype, color = values)) + geom_point(size = 4) + guides(color = FALSE) + ylab(element_blank()) + xlab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

drift.index.Subtypes.ast <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers.ast, species.1 = gene.lists.pos$surface.markers, species.2 = gene.lists.pos$cave.markers)
DI.ast <- tibble(Subtype = c(names(unlist(drift.index.Subtypes.ast))), values = c(unlist(drift.index.Subtypes.ast)))
# ggplot(DI.ast, aes(y = Subtype, x = values, label = Subtype, color = values)) + geom_point(size = 4) + guides(color = FALSE) + ylab(element_blank()) + xlab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


## Plot Subclusters
drift.index.SubclusterType <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub)
DI.sub <- tibble(SubclusterType = c(names(unlist(drift.index.SubclusterType))), values = c(unlist(drift.index.SubclusterType)))
# ggplot(DI.sub, aes(x = SubclusterType, y = values, label = SubclusterType, color = values)) + geom_point(size = 4) + guides(color = FALSE) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

drift.index.SubclusterTypes.ast <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers.ast.sub, species.1 = gene.lists.pos$surface.markers.sub, species.2 = gene.lists.pos$cave.markers.sub)
DI.ast.sub <- tibble(SubclusterType = c(names(unlist(drift.index.SubclusterTypes.ast))), values = c(unlist(drift.index.SubclusterTypes.ast)))
# ggplot(DI.ast.sub, aes(y = SubclusterType, x = values, label = Subtype, color = values)) + geom_point(size = 4) + guides(color = FALSE) + ylab(element_blank()) + xlab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
 
# OK, make a df with the Subtype as a column
# Plot subtype values ontop!

index <- unique(tibble(Subtype = hypo.integrated@meta.data$integrated_Subtype, SubclusterType = hypo.integrated@meta.data$integrated_SubclusterType))
DI.sub$Subtype <- index$Subtype[match(DI.sub$SubclusterType, index$SubclusterType)]
DI.sub$Subtype <- factor(DI.sub$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

index.ast <- unique(tibble(Subtype = hypo.integrated.ast@meta.data$integrated_Subtype, SubclusterType = hypo.integrated.ast@meta.data$integrated_SubclusterType))
DI.ast.sub$Subtype <- index.ast$Subtype[match(DI.ast.sub$SubclusterType, index.ast$SubclusterType)]
DI.ast.sub$Subtype <- factor(DI.ast.sub$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

DI.plot <- ggplot(DI.sub, aes(x = Subtype, y = values, colour = Subtype, fill = Subtype, group = Subtype)) + geom_boxplot(width = 0.5) + geom_point(size = 1, color = "black") + guides(color = FALSE, fill = FALSE) + xlab(element_blank()) + ylab("Drift Index") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DI.plot <- DI.plot + geom_point(data = DI, aes(x = Subtype, y = values), size = 3, shape = 18, color = "red", position = position_nudge(x=0.4)) + theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 8), axis.text.x = element_text(size = 8), axis.title = element_text(size = 10))

png("Figures/Drift_Index_species_subtype.png", width = 8, height = 5, units = "in", res = 250)
DI.plot
dev.off()

DI.plot.ast <- ggplot(data = DI.ast.sub, aes(x = Subtype, y = values, colour = Subtype, fill = Subtype, group = Subtype)) + geom_boxplot(width = 0.5) + geom_point(size = 1, color = "black") + guides(color = FALSE, fill = FALSE) + xlab(element_blank()) + ylab("Drift Index") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DI.plot.ast <- DI.plot.ast + geom_point(data = DI.ast, aes(x = Subtype, y = values), size = 3, shape = 18, color = "blue", position = position_nudge(x=0.4)) + theme(axis.text.y = element_text(hjust = 1, vjust = 0.5, size = 8), axis.text.x = element_text(size = 8), axis.title = element_text(size = 10))

png("Figures/Drift_Index_morphs_subtype.png", width = 10, height = 6, units = "in", res = 250)
DI.plot.ast
dev.off()


# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())][4:9]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files()[grep("GO", list.files())][4:9]
names(go_lists)
go_lists[[7]] <- Reduce(union, list(go_lists[[3]], go_lists[[4]], go_lists[[5]]))


## Plot Subtypes

test <- lapply(go_lists, function(x) c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers, species.1 = gene.lists.pos$zebrafish.markers, species.2 = gene.lists.pos$astyanax.markers, subset = x, invert = F))))
names(test) <- names(go_lists)

test <- Reduce(cbind, test)
colnames(test) <- names(go_lists)

DI2 <- tibble(cell_type = c(names(unlist(drift.index.SubclusterType))), all = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, invert = F))), TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = go_lists[[2]], invert = F))), NP_NTS = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = go_lists[[7]], invert = F))))

DI2 <- DI2[!is.na(DI2$NP_NTS),]

DI3 <- melt(DI2)

ggplot(DI3, aes(x = variable, y = value, group = variable, color = variable)) + geom_jitter(size = 1) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_manual(values = c("blue", "red", "green")) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means()


## Plot Subclusters

test <- lapply(go_lists, function(x) c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = x, invert = F))))

test <- data.frame(axon = test[[2]], TF = test[[3]], gaba = test[[4]], glut = test[[5]], NP = test[[6]])
test$cell_type <- row.names(test)
test <- melt(test)

DI2 <- tibble(cell_type = Reduce(intersect, list(names(gene.lists.pos$zebrafish.markers.sub), names(gene.lists.pos$astyanax.markers.sub), names(gene.lists.pos$conserved.markers.sub))), TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = go_lists[[3]], invert = F))), NPs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = go_lists[[6]], invert = F))), Non_TF_NPs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = c(go_lists[[3]], go_lists[[6]]), invert = T))))

DI2 <- DI2[!is.na(DI2$NPs),]


DI3 <- melt(DI2)

DI3$cell_type <- factor(DI3$cell_type, levels = levels(sort(DI3$cell_type)))

ggplot(DI3, aes(x = cell_type, y = value, group = variable, color = variable)) + geom_point(size = 2) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + geom_line()

DI4 <- DI3[grep("GABA_1", DI3$cell_type),]
DI4$variable <- factor(DI4$variable, levels = c("all", "NP_NTS", "TFs"))
DI4$cell_type <- factor(DI4$cell_type, levels = c("GABA_1_0", "GABA_1_1", "GABA_1_2", "GABA_1_3", "GABA_1_4", "GABA_1_5", "GABA_1_6", "GABA_1_7", "GABA_1_8", "GABA_1_9", "GABA_1_10", "GABA_1_11", "GABA_1_12", "GABA_1_13", "GABA_1_14"))

ggplot(DI4, aes(x = cell_type, y = value, group = variable, color = variable)) + geom_point(size = 2) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), panel.background = element_blank(), axis.line = element_line(color = "black")) + geom_line() + scale_color_viridis_d()

ggplot(DI3, aes(x = variable, y = value, group = variable, color = variable)) + geom_boxplot(outlier.color = NA) + geom_jitter(size = 1) + xlab(element_blank()) + ylab("Cellular Drift Index") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means(method = "wilcox")




