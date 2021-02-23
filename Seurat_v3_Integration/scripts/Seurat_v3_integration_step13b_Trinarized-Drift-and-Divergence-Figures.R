library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(viridis)
library(ggpubr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

a = 1.5
b = 2
f = 0.1

## Load trinarized gene lists

trinarized.genes <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,".rds", sep = ""))
trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# Add species-specific subclusters
norm.cluster.filtered.zeb <- lapply(c(trinarized.exp$subcluster.zebrafish, trinarized.exp$specific.zebrafish), function(x) names(x))
norm.cluster.filtered.ast <- lapply(c(trinarized.exp$subcluster.astyanax, trinarized.exp$specific.astyanax), function(x) names(x))

## Divergence (dT) plots
## Figure 2d-e, S4c

dT.list <- readRDS(file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

branch.levels <- levels(dT.list[[2]]$branch)
branch.levels <- c(branch.levels, "Danio rerio")

dT.list[[1]]$branch <- factor(dT.list[[1]]$branch, levels = branch.levels)
dT.list[[2]]$branch <- factor(dT.list[[2]]$branch, levels = branch.levels)

plot.zeb <- ggplot(dT.list[[1]], aes(x = dT, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour") + theme_classic()
plot.ast <- ggplot(dT.list[[2]], aes(x = dT, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour") + theme_classic()

plot.zeb <- plot.zeb + guides(color = guide_legend(title = "LCA of paralogs"))
plot.ast <- plot.ast + guides(color = guide_legend(title = "LCA of paralogs"))


## Compare gene pairs which are common to zebrafish and cavefish
## Compare dTs for cavefish and zebrafish, and correlate binarized gene expression patterns

library(plyr)

dTzeb.2 <- plyr::match_df(dT.list[[1]], dT.list[[2]], on = c("gene1", "gene2"))
dTast.2 <- plyr::match_df(dT.list[[2]], dT.list[[1]], on = c("gene1", "gene2"))
dT.combined <- merge(dTast.2, dTzeb.2, by = c("gene1", "gene2")) # all = T for merging everything and making up rows, below still works
colnames(dT.combined) <- c("gene1", "gene2", "dT.ast", "branch.ast", "dT.zeb", "branch.zeb")
# x is ast, y is zeb
all.genes <- ggplot(dT.combined, aes(x = dT.x, y = dT.y, color = branch.x, label = paste(gene1, gene2, sep = "-"))) + geom_point(colour = "black", size = 0.5) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + theme_classic() + ylab("dT zebrafish") + xlab("dT Mexican tetra") + ggtitle("All gene pairs", subtitle = "Pearson correlation = 0.54")
otophysi <- ggplot(dT.combined[dT.combined$branch.x == "Otophysi" | dT.combined$branch.y == "Otophysi",], aes(x = dT.x, y = dT.y)) + geom_point(colour = "#FF61C9", size = 0.5) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + theme_classic() + ylab("dT zebrafish") + xlab("dT Mexican tetra") + ggtitle("LCA Otophysi gene pairs", subtitle = "Pearson correlation = 0.45")

all.genes + otophysi




## Make Paralog divergence supplemental figure S7

##### Correlate binarized cell types between species for each pair

# Create matrices of normalized expression for each gene across integrated_SubclusterTypes

zeb.int <- Reduce(cbind, trinarized.genes$subcluster.zebrafish)
colnames(zeb.int) <- names(trinarized.genes$subcluster.zebrafish)
# row.names(zeb.int) <- row.names(trinarized.exp$subcluster.zebrafish[[1]])

ast.int <- Reduce(cbind, trinarized.genes$subcluster.astyanax)
colnames(ast.int) <- names(trinarized.genes$subcluster.astyanax)
# row.names(ast.int) <- row.names(trinarized.exp$subcluster.astyanax[[1]])

dTcombined <- lapply(unique(dT.combined$branch.x), function(x) dT.combined[dT.combined$branch.x == x,])
names(dTcombined) <- unique(dT.combined$branch.x)

# Match the gene pairs from the dTcombined lists from each .int matrix and cor() them

correlations <- lapply(dTcombined, function(x) apply(x, 1, function(x) cor(unlist(c(zeb.int[as.character(x[1]),intersect(colnames(zeb.int), colnames(ast.int))], zeb.int[as.character(x[2]),intersect(colnames(zeb.int), colnames(ast.int))])), unlist(c(ast.int[as.character(x[1]),intersect(colnames(zeb.int), colnames(ast.int))], ast.int[as.character(x[2]),intersect(colnames(zeb.int), colnames(ast.int))])))))

lapply(correlations, function(x) mean(x))

correlations2 <- data.frame(branch = c(rep(names(correlations)[[1]], length(correlations[[1]])), rep(names(correlations)[[2]], length(correlations[[2]])), rep(names(correlations)[[3]], length(correlations[[3]])), rep(names(correlations)[[4]], length(correlations[[4]])), rep(names(correlations)[[5]], length(correlations[[5]])), rep(names(correlations)[[6]], length(correlations[[6]])), rep(names(correlations)[[7]], length(correlations[[7]])), rep(names(correlations)[[8]], length(correlations[[8]]))), value = unlist(correlations))

correlations2$branch <- factor(correlations2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))

anno_df <- compare_means(value ~ branch, data = correlations2)

# Figure S7a
corr.plot <- ggplot(correlations2, aes(x = branch, y = value)) + geom_jitter(color = "grey65", shape = 21) + geom_boxplot(fill = "transparent", outlier.color = NA) + theme_classic() + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + stat_compare_means()


## Do it for binarized expression (there or not)

zeb.int.b <- ifelse(zeb.int > .95, 1, 0)
ast.int.b <- ifelse(ast.int > .95, 1, 0)

correlations.b <- lapply(dTcombined, function(x) apply(x, 1, function(x) cor(unlist(c(zeb.int.b[as.character(x[1]),intersect(colnames(zeb.int.b), colnames(ast.int.b))], zeb.int.b[as.character(x[2]),intersect(colnames(zeb.int.b), colnames(ast.int.b))])), unlist(c(ast.int.b[as.character(x[1]),intersect(colnames(zeb.int.b), colnames(ast.int.b))], ast.int.b[as.character(x[2]),intersect(colnames(zeb.int.b), colnames(ast.int.b))])))))

lapply(correlations.b , function(x) mean(x))

correlations.b.2 <- data.frame(branch = c(rep(names(correlations)[[1]], length(correlations.b[[1]])), rep(names(correlations)[[2]], length(correlations.b[[2]])), rep(names(correlations)[[3]], length(correlations.b[[3]])), rep(names(correlations)[[4]], length(correlations.b[[4]])), rep(names(correlations)[[5]], length(correlations.b[[5]])), rep(names(correlations)[[6]], length(correlations.b[[6]])), rep(names(correlations)[[7]], length(correlations.b[[7]])), rep(names(correlations)[[8]], length(correlations.b[[8]]))), value = unlist(correlations.b))

correlations.b.2$branch <- factor(correlations.b.2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))

anno_df <- compare_means(value ~ branch, data = correlations.b.2)

# Figure S7b
s7b <- ggplot(correlations.b.2[correlations.b.2$branch != "Opisthokonta",], aes(x = branch, y = value)) + geom_jitter(color = "grey65", shape = 21) + geom_boxplot(fill = "transparent", outlier.color = NA) + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + stat_compare_means()

s7b <- s7b + theme_classic() + stat_compare_means()



## Plot dT.ancestral and NF scores

dT.NF <- dT.list[[4]]

## Make some plots

# Again, a high dT or ancestral dT in one species does not predict a high dT in the other
ggplot(dT.NF, aes(x = dT.zeb, y = dT.ast)) + geom_point() + theme_classic()
ggplot(dT.NF, aes(x = dT.zeb.ancestral, y = dT.ast.ancestral)) + geom_point() + theme_classic()

# NF scores for gene pairs are poorly correlated, but most likely due to the nature of how ancestral cell types are calculated
ggplot(dT.NF, aes(x = NF.gene1.zeb, y = NF.gene2.zeb)) + geom_point() + theme_classic()
ggplot(dT.NF, aes(x = NF.gene1.ast, y = NF.gene2.ast)) + geom_point() + theme_classic()

ggplot(dT.NF, aes(x = NF.gene1.zeb+NF.gene2.zeb, y = NF.gene1.ast+NF.gene2.ast)) + geom_point() + theme_classic()

# No relationship between dT and combined NF scores
ggplot(dT.NF, aes(x = dT.zeb.ancestral, y = NF.gene1.zeb+NF.gene2.zeb)) + geom_point() + theme_classic()
ggplot(dT.NF, aes(x = dT.ast.ancestral, y = NF.gene1.ast+NF.gene2.ast)) + geom_point() + theme_classic()


ggplot(dT.NF, aes(x = branch, y = NF.gene1.zeb+NF.gene2.zeb)) + geom_jitter() + geom_boxplot() + theme_classic()
ggplot(dT.NF, aes(x = branch, y = NF.gene1.ast+NF.gene2.ast)) + geom_jitter() + geom_boxplot() + theme_classic()

ggplot(dT.NF, aes(colour = branch, x = NF.gene1.zeb+NF.gene2.zeb)) + stat_ecdf() + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour") + theme_classic()
ggplot(dT.NF, aes(colour = branch, x = NF.gene1.ast+NF.gene2.ast)) + stat_ecdf() + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour") + theme_classic()


# Melt the NF/dT values into single column for easier plotting
ggplot(melt(dT.NF[,9:13]), aes(x = branch, y = value)) + geom_jitter() + geom_boxplot() + theme_classic()
ggplot(melt(dT.NF[,9:13]), aes(x = value, color = branch)) + stat_ecdf(size = 1, geom = "step") + theme_classic() + scale_colour_viridis(discrete = T)

dT.plot <- ggplot(melt(dT.NF[,c(7:8,13)]), aes(x = value, color = branch)) + stat_ecdf(size = 1, geom = "step") + theme_classic() + scale_colour_viridis(discrete = T)
dT.ancestral.plot <- ggplot(melt(dT.NF[,c(3,5,13)]), aes(x = value, color = branch)) + stat_ecdf(size = 1, geom = "step") + theme_classic() + scale_colour_viridis(discrete = T)

dT.plot + dT.ancestral.plot




ggplot(dT.NF, aes(x = branch, y = dT.ast.ancestral)) + geom_jitter() + geom_boxplot() + theme_classic()
ggplot(dT.NF, aes(x = branch, y = dT.zeb.ancestral)) + geom_jitter() + geom_boxplot() + theme_classic()


ggplot(dT.NF, aes(x = dT.zeb, y = NF.gene1.zeb+NF.gene2.zeb, colour = branch)) + geom_point() + theme_classic()
ggplot(dT.NF, aes(x = dT.ast, y = NF.gene1.ast+NF.gene2.ast, colour = branch)) + geom_point() + theme_classic()


ggplot(dT.NF, aes(y = NF.gene1.zeb/NF.gene2.zeb, x = branch, colour = branch)) + geom_jitter() + geom_boxplot() + theme_classic()

ggplot(dT.NF, aes(y = NF.gene1.ast, x = NF.gene2.ast, colour = branch)) + geom_point() + theme_classic()


