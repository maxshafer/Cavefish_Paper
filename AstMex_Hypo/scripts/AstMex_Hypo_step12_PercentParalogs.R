library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(data.table)
library(fdrtool)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/")

## Take sets of conserved, and species specific marker genes (start with Subtypes)
## For each species, take each gene (apply), and query mart database for paralogs
## Ask if any of the paralogs are in species.2 list
## Append species.1 list with paralog expressed (if more then 1??), plus LCA


### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Otophysi_Paralogs/mart_export_Danio.txt", sep = "\t", head = TRUE)
mart[[2]] <- read.csv("~/Downloads/mart_export_Astyanax_102.txt", sep = "\t", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")


## Load marker gene lists

gene.lists <- readRDS(file = "/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/drift_gene_lists.rds")



## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4)] <- lapply(gene.lists[c(1,4)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$zebrafish_avg_logFC > 0 & x[[y]]$astyanax_avg_logFC > 0,]))
gene.lists.pos[c(2,3,5,6,8,9,11,12,14,15,17,18)] <- lapply(gene.lists[c(2,3,5,6,8,9,11,12,14,15,17,18)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))
gene.lists.pos[c(7,10,13,16)] <- lapply(gene.lists[c(7,10,13,16)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

### Association between whether a gene is a specieis specific marker, and it being a paralog of a marker gene in either species, or conserved between species.

## Identify all paralogs of all marker genes union(paralogs of conserved, paralogs of species.1, paraglos of species.2)
## Ask what % of species specific marker genes are in the above list

calcParalog <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, mart.1 = mart[[1]], mart.2 = mart[[2]], ngenes.1 = 25271, ngenes.2 = 25271, i = i) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	
	paralog.con <- unique(mart.2$Cave.fish.paralogue.associated.gene.name[match(row.names(conserved[[i]]), mart.2$Gene.name)])
	
	paralog.1 <- unique(mart.2$Cave.fish.paralogue.associated.gene.name[match(setdiff(row.names(species.1[[i]]), row.names(conserved[[i]])), mart.2$Gene.name)])
	
	paralog.2 <- unique(mart.2$Cave.fish.paralogue.associated.gene.name[match(setdiff(row.names(species.2[[i]]), row.names(conserved[[i]])), mart.2$Gene.name)])
	paralog.union <- union(paralog.con, union(paralog.1, paralog.2))
	
	genes.1 <- setdiff(row.names(species.1[[i]]), row.names(conserved[[i]]))
	genes.2 <- setdiff(row.names(species.2[[i]]), row.names(conserved[[i]]))
	
	a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.2)])
	b1 <- length(genes.1) - a1
	c1 <- length(union(paralog.con, paralog.2)) - a1
	d1 <- ngenes.1 - b1

	a2 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.1)])
	b2 <- length(genes.1) - a2
	c2 <- length(union(paralog.con, paralog.1)) - a2
	d2 <- ngenes.2 - b2
	
	a3 <- length(row.names(conserved[[i]]))
	b3 <- length(genes.1)
	c3 <- length(genes.2)
	d3 <- ngenes.1 - sum(b3, c3)
	
	vec <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2, a3=a3, b3=b3, c3=c3, d3=d3)
	return(vec)
}



paralog.numbers <- lapply(seq_along(gene.lists.pos[[13]]), function(x) calcParalog(conserved = gene.lists.pos[[13]], species.1 = gene.lists.pos[[14]], species.2 = gene.lists.pos[[15]], i = x))


fisher.results.1 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[[13]])
fisher.results.2 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(gene.lists.pos[[13]])
fisher.results.3 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2)))
fisher.results.3 <- do.call(rbind, fisher.results.3)
row.names(fisher.results.3) <- names(gene.lists.pos[[13]])

subtypes <- data.frame(do.call(rbind, paralog.numbers))
row.names(subtypes) <- names(gene.lists.pos[[13]])
subtypes$percent.para.1 <- unlist(apply(subtypes, 1, function(x) x[1]/x[2]))*100
subtypes$percent.para.2 <- unlist(apply(subtypes, 1, function(x) x[5]/x[6]))*100
# subtypes$fisher.pval <- as.numeric(fisher.results[,1])
subtypes$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
subtypes$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
subtypes$fisher.pval.1 <- as.numeric(fisher.results.1[,1])
subtypes$fisher.pval.2 <- as.numeric(fisher.results.2[,1])
subtypes$Subtypes <- row.names(subtypes)
subtypes$driftindex <- DI$values

subtypes$Subtypes <- factor(subtypes$Subtypes, levels = c("Endothelial", "Erythrocytes", "Ciliated/Ventricular", "Ependymal", "Progenitors_1", "Progenitors_2", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes_1", "Oligodendrocytes_2", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "GABA_6", "GABA_7", "GABA_Prdx1_1", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Galanin", "Otpa/b_1", "Otpa/b_2", "Otpa/b_3", "Lymphatic", "Tcells", "Bcells", "Mast_cells", "Thrombocytes", "Neutrophils", "Macrophages", "Microglia"))

subtypes.2 <- melt(subtypes[,13:18])

#### For subclusters!

paralog.numbers.sub <- lapply(seq_along(gene.lists.pos[[16]]), function(x) calcParalog(conserved = gene.lists.pos[[16]], species.1 = gene.lists.pos[[17]], species.2 = gene.lists.pos[[18]], i = x))

fisher.results.1 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[[16]])
fisher.results.2 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(gene.lists.pos[[16]])

subclusters <- data.frame(do.call(rbind, paralog.numbers.sub))
row.names(subclusters) <- names(gene.lists.pos[[16]])
subclusters$percent.para.1 <- unlist(apply(subclusters, 1, function(x) x[1]/x[2]))*100
subclusters$percent.para.2 <- unlist(apply(subclusters, 1, function(x) x[5]/x[6]))*100
# subclusters$fisher.pval <- as.numeric(fisher.results[,1])
subclusters$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
subclusters$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
subclusters$fisher.pval.1 <- as.numeric(fisher.results.1[,1])
subclusters$fisher.pval.2 <- as.numeric(fisher.results.2[,1])
subclusters$Subclusters <- row.names(subclusters)
subclusters$driftindex <- DI.sub$values
subclusters$Subtype <- hypo.integrated.ast@meta.data$Subtype[match(subclusters$Subclusters, hypo.integrated.ast@meta.data$SubclusterType)]

subclusters$Subtype <- factor(subclusters$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated/Ventricular", "Ependymal", "Progenitors_1", "Progenitors_2", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes_1", "Oligodendrocytes_2", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "GABA_6", "GABA_7", "GABA_Prdx1_1", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Galanin", "Otpa/b_1", "Otpa/b_2", "Otpa/b_3", "Lymphatic", "Tcells", "Bcells", "Mast_cells", "Thrombocytes", "Neutrophils", "Macrophages", "Microglia"))

subclusters.2 <- melt(subclusters[,13:19])


# ## OK have several dataframes:
# ## subtypes, subtypes.2(melted), subclusters, subclusters.2(metled)

# subtypes.3 <- data.frame(Subtypes = rep(subtypes$Subtypes, 2), driftindex = rep(subtypes$driftindex, 2), species = c(rep("drerio", 25), rep("amexicanus", 25)), percent.para = c(subtypes$percent.para.1, subtypes$percent.para.2), fisher.odds = c(subtypes$fisher.odds.1, subtypes$fisher.odds.2))


### PLOTS

## Make DI by Paralog (sum?? or mean??)

m <- lm((percent.para.2 + percent.para.1)/2 ~ driftindex, subclusters)
r2 <- summary(m)$r.squared

r <- cor((subclusters$percent.para.2 + subclusters$percent.para.1)/2, subclusters$driftindex)

para.drift <- ggplot(subclusters, aes(x = (percent.para.2 + percent.para.1)/2, y = driftindex, color = Subtype, label = Subclusters)) + geom_point() + guides(color = F) + geom_smooth(method = lm, se = T, color = "black", fill = "grey85") + ylab("Drift Index") + xlab("Mean of % Paralogs per species")

para.drift <- para.drift + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + geom_text(x = 2.5, y = 0.6, label = r, parse = T, color = "black")


## Make Lollipop plot for expected vs observed enrichment

percent.para <- ggplot(subtypes) + geom_segment( aes(x = percent.para.1/fisher.odds.1, xend = percent.para.1, y = Subtypes, yend = Subtypes), colour = "black", position = position_nudge(y=0.2)) + geom_point( aes(x = percent.para.1/fisher.odds.1, y = Subtypes, colour = Subtypes), alpha = 0.5, size = 2, position = position_nudge(y=0.2)) + geom_point( aes(x = percent.para.1, y = Subtypes, colour = Subtypes), size = 2, position = position_nudge(y=0.2)) + geom_segment( aes(x = percent.para.2/fisher.odds.2, xend = percent.para.2, y = Subtypes, yend = Subtypes), colour = "red", position = position_nudge(y=-0.2)) + geom_point( aes(x = percent.para.2/fisher.odds.2, y = Subtypes, colour = Subtypes), alpha = 0.5, size = 2, position = position_nudge(y=-0.2)) + geom_point( aes(x = percent.para.2, y = Subtypes, colour = Subtypes), size = 2, position = position_nudge(y=-0.2)) + guides(colour = F) + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + ylab("") + xlab("Expected vs \n Observed % Paralogs") + coord_flip()




# subtypes <- melt(subtypes[,9:12])

ggplot(subtypes, aes(x = percent.para.2 + percent.para.1, y = driftindex, color = Subtypes, label = Subtypes)) + geom_point() + geom_text() + guides(color = FALSE) + geom_smooth(method = lm, se = T, color = "black", fill = "grey85")
ggplot(subtypes, aes(x = percent.para.1, y = percent.para.2, color = Subtypes, label = Subtypes)) + geom_point() + geom_text() + guides(color = FALSE)



ggplot(subtypes.2[subtypes.2$variable == "percent.para.1" | subtypes.2$variable == "percent.para.2",], aes(Subtypes,x = variable, y = value, color = Subtypes, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F) + ylab("Percentage Paralog") + xlab("Species") + scale_x_discrete(labels = c("D. rerio", "A. mexicanus")) + theme(axis.text.x = element_text(face = "italic"))

ggplot(subtypes.2[subtypes.2$variable == "fisher.odds.1" | subtypes.2$variable == "fisher.odds.2",], aes(Subtypes,x = variable, y = value, color = Subtypes, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F)
ggplot(subtypes.2[subtypes.2$variable == "driftindex",], aes(Subtypes,x = variable, y = value, color = Subtypes, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F)


fisher.test(matrix(data = apply(do.call(rbind, paralog.numbers), 2, function(x) sum(x))[5:8], ncol = 2, nrow = 2))



## FOR subclusters

ggplot(subclusters, aes(x = percent.para.2 + percent.para.1, y = driftindex, color = Subclusters, label = Subclusters)) + geom_point() + guides(color = FALSE) + geom_smooth(method = lm, se = T, color = "black", fill = "grey85")
ggplot(subclusters, aes(x = percent.para.1, y = driftindex, color = Subclusters, label = Subclusters)) + geom_point() + guides(color = FALSE)


ggplot(subclusters.2[subclusters.2$variable == "percent.para.1" | subclusters.2$variable == "percent.para.2",], aes(Subclusters,x = variable, y = value, color = Subclusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 1, position = position_jitter(0.1)) + guides(color = F)

ggplot(subclusters.2[subclusters.2$variable == "fisher.odds.1" | subclusters.2$variable == "fisher.odds.2",], aes(Subclusters, x = variable, y = value, color = Subclusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 1, position = position_jitter(0.1)) + guides(color = F)

ggplot(subclusters.2[subclusters.2$variable == "driftindex",], aes(Subclusters,x = variable, y = value, color = Subclusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 1, position = position_jitter(0.1)) + guides(color = F)


fisher.test(matrix(data = apply(do.call(rbind, paralog.numbers.sub), 2, function(x) sum(x))[5:8], ncol = 2, nrow = 2))



 
## save plots

percent.para <- ggplot(subtypes.2[subtypes.2$variable == "percent.para.1" | subtypes.2$variable == "percent.para.2",], aes(Subtypes,x = variable, y = value, color = Subtypes, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F) + ylab("Percentage Paralog") + xlab("Species") + scale_x_discrete(labels = c("D. rerio", "A. mexicanus")) + theme(axis.text.x = element_text(face = "italic"))

fisher.odds <- ggplot(subtypes.2[subtypes.2$variable == "fisher.odds.1" | subtypes.2$variable == "fisher.odds.2",], aes(Subtypes,x = variable, y = value, color = Subtypes, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F) + ylab("Fisher Odds") + xlab("Species") + scale_x_discrete(labels = c("D. rerio", "A. mexicanus")) + theme(axis.text.x = element_text(face = "italic"))





## Make Patchwork plot for Figure 3
library(patchwork)

# Takes plots from the PercentParalogs.R, Drift_Index.R, and expression_divergence_v2.R scripts

( (DI.plot | para.drift) + plot_layout(ncol = 2, widths = c(3,1)) ) / ( (percent.para | plot.id) + plot_layout(ncol = 2, widths = c(3,1)) )



