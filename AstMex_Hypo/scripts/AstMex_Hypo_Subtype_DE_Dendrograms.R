library(Seurat)
library(tidyr)
library(dendextend)
library(dplyr)
library(ape)
library(circlize)
library(ggplot2)
library(scales)
library(tibble)


# Load seurat object (astmex)
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/")
load("AstMex_64k.Robj")

###################################################
# Generate proportional stacked bar plot
###################################################

# Make colour palettes

ids <- c("Subtype", "SubclusterType")

my_colour_palette <- list()
for (i in 1:length(ids)) {
	hypo <- SetAllIdent(hypo, id = ids[[i]])
	colours <- hue_pal()(length(levels(hypo@ident)))
	names(colours) <- levels(hypo@ident)
	colours <- colours[dendrograms[[i]] %>% labels]
	my_colour_palette[[i]] <- colours
}

# Make tables of cell type proportions
prop.table <- list()
prop.table[[1]] <- table(hypo@meta.data$Subtype, hypo@meta.data$species_morph)
prop.table[[2]] <- table(hypo@meta.data$SubclusterType, hypo@meta.data$species_morph)
prop.table[[3]] <- table(hypo@meta.data$Subtype, hypo@meta.data$species)

prop.table <- lapply(prop.table, function(x) as.data.frame(t(apply(x, 1, function(y) {y/sum(y)}))))

for (i in 1:length(prop.table)) {
	prop.table[[i]]$cell_type <- row.names(prop.table[[i]])
	prop.table[[i]] <- prop.table[[i]][,c(5,1,2,3,4)]
	prop.table[[i]] <- prop.table[[i]] %>% gather(species_morph, freq, choy_surface:molino_cave:pachon_cave:tinaja_cave)
	prop.table[[i]]$cell_type <- as.factor(prop.table[[i]]$cell_type)
}

prop.plots <- list()
for (i in 1:length(prop.table)) {
	prop.plots[[i]] <- ggplot(prop.table[[i]], aes(x=cell_type, y=freq, fill=species_morph)) + geom_bar(stat="identity") + guides(fill = FALSE) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_text(size = 8), axis.title.x = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = c("springgreen4", "goldenrod1", "lightgoldenrod1", "darkorange1")) + ylab("% species \n morph")
}

###################################################
# Generate diff and snp gene plots
###################################################

#

species.markers <- readRDS("AstMex_Hypo_markers.species.Subtype.list.rds")

# Diff genes plot

species.markers.pval <- lapply(species.markers, function(x) x[x$p_val_adj < 0.05,])

subtype.diff.counts <- unlist(lapply(species.markers.pval, function(x) nrow(x)))

subtype.diff.counts <- as.data.frame(subtype.diff.counts)
subtype.diff.counts$cell_type <- rownames(subtype.diff.counts)
colnames(subtype.diff.counts) <- c("diff.counts", "cell_type")

subtype.diff.counts <- subtype.diff.counts[subtype.diff.counts$cell_type %in% levels(hypo@meta.data$Subtype),]

diff.plot <- ggplot(subtype.diff.counts, aes(x=cell_type, y = diff.counts)) + geom_bar(stat="identity", fill = "red") + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 8), plot.background = element_rect(fill = "transparent", color = NA)) + ylab("# of \n DE genes")


###################################################
# Generate cluster similarity plot
###################################################

# Generate phylogenetic distances (similarity between clusters)

dendrograms <- hypo@cluster.tree

phylo <- list()
phylo <- lapply(dendrograms, function(x) cophenetic.phylo(as.phylo(x)))
phylo <- as.data.frame(cophenetic.phylo(as.phylo(dendrograms[[3]])))

phylo <- as.matrix(phylo)

phylo.dist <- diag(phylo[37:72,1:36])
names(phylo.dist) <- row.names(phylo[37:72,])
names(phylo.dist) <- sub('.*e ', "", names(phylo.dist))
phylo.dist <- data.frame(cell_type = names(phylo.dist), dist = phylo.dist)

phylo.plot <- ggplot(phylo.dist, aes(y = dist, x = cell_type)) + geom_bar(stat = "identity", fill = "black") + guides(fill = FALSE) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_text(size = 8), plot.background = element_rect(fill = "transparent", color = NA)) + ylab("Dissimilarity")

###################################################
# Generate combine plots!
###################################################

plots <- phylo.plot / diff.plot / prop.plots[[1]]
plots & theme(axis.text = element_text(size = 8))

plot_grid(phylo.plot, diff.plot, prop.plots[[1]], nrow = 3, rel_heights = c(.5, .5, 1.5), align = "v")
plot_grid(phylo.plot, diff.plot, nrow = 2, rel_heights = c(1, 1), align = "v")

ggsave("Figures/AstMex_barplots_subtype.png", units = "in", height = 7, width = 12)





