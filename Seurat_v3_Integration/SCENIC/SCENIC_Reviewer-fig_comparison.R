library(reshape2)
library(Seurat)
library(viridis)
library(ggrepel)
library(scales)
library(ggsignif)
library(ggpubr)
library(patchwork)


# Load old GENIE3 files from all the comps, for ast and zeb

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC copy")

gene.lists <- c("nps", "nts", "synaptic", "ion")

zeb.lists <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
ast.lists <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
zeb.modules <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
ast.modules <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))

names(zeb.lists) <- gene.lists
names(ast.lists) <- gene.lists
names(zeb.modules) <- gene.lists
names(ast.modules) <- gene.lists

## Make functions for extracting/plotting data
## load functions
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/SCENIC_functions.R")

## Load GO_lists and make plot data and plots

GO_lists <- lapply(gene.lists, function(x) union(unique(zeb.lists[[x]]$Target), unique(ast.lists[[x]]$Target)))
names(GO_lists) <- gene.lists
GO_data <- lapply(seq_along(GO_lists), function(y) lapply(seq_along(GO_lists[[y]]), function(x) extractPlotData(zeb = zeb.lists[[y]], ast = ast.lists[[y]], cor.zeb = zeb.modules[[y]], cor.ast = ast.modules[[y]], np = GO_lists[[y]][x])))

names(GO_data) <- names(GO_lists)
for (i in 1:length(GO_data)) {
  names(GO_data[[i]]) <- GO_lists[[i]]
}

GO_data.old <- GO_data

# Load new GENIE3 files from all the comps, for ast and zeb

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC")

# Load GENIE3 files from all the comps, for ast and zeb

gene.lists <- c("nps", "nts", "synaptic", "ion")

zeb.lists <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
ast.lists <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
zeb.modules <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
ast.modules <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))

names(zeb.lists) <- gene.lists
names(ast.lists) <- gene.lists
names(zeb.modules) <- gene.lists
names(ast.modules) <- gene.lists

## Make functions for extracting/plotting data
## load functions
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/SCENIC_functions.R")

## Load GO_lists and make plot data and plots

GO_lists <- lapply(gene.lists, function(x) union(unique(zeb.lists[[x]]$Target), unique(ast.lists[[x]]$Target)))
names(GO_lists) <- gene.lists
GO_data <- lapply(seq_along(GO_lists), function(y) lapply(seq_along(GO_lists[[y]]), function(x) extractPlotData(zeb = zeb.lists[[y]], ast = ast.lists[[y]], cor.zeb = zeb.modules[[y]], cor.ast = ast.modules[[y]], np = GO_lists[[y]][x])))

names(GO_data) <- names(GO_lists)
for (i in 1:length(GO_data)) {
  names(GO_data[[i]]) <- GO_lists[[i]]
}

GO_data.new <- GO_data


### Compare the results

# Each GO_data is a list of a list of data.frames
# Need to only compare TFs that are present in both datasets (then can also compare the correlations between species)
# need to find the intersect of TFs, then only compare the weights (for zeb and ast) for those TFs between old and new
TF.old <- unique(Reduce(rbind, lapply(GO_data.old, function(x) Reduce(rbind, x)))$TF)
TF.new <- unique(Reduce(rbind, lapply(GO_data.new, function(x) Reduce(rbind, x)))$TF)

intersect <- intersect(TF.old, TF.new)

# hmmm, even with the intersect, there are different TFs in each data.frame - need to subset by the intersection of the TFs
cor(GO_data.old$nps$oxt$weight_zeb[GO_data.old$nps$oxt$TF %in% intersect], GO_data.new$nps$oxt$weight_zeb[GO_data.new$nps$oxt$TF %in% intersect])

cor(GO_data.old$nps$oxt$weight_ast[GO_data.old$nps$oxt$TF %in% intersect], GO_data.new$nps$oxt$weight_ast[GO_data.new$nps$oxt$TF %in% intersect])

cor <- list()
cor[[1]] <- list()
cor[[2]] <- list()
names(cor) <- c("zebrafish_grn", "astyanax_grn")
nps.intersect <- intersect(names(GO_data.old$nps), names(GO_data.new$nps))
for (i in 1:length(nps.intersect)) {
  tf.intersect <- intersect(GO_data.old$nps[[nps.intersect[[i]]]]$TF, GO_data.new$nps[[nps.intersect[[i]]]]$TF)
  cor[[1]][[i]] <- cor(GO_data.old$nps[[nps.intersect[[i]]]]$weight_zeb[GO_data.old$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect], GO_data.new$nps[[nps.intersect[[i]]]]$weight_zeb[GO_data.new$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect])
  cor[[2]][[i]] <- cor(GO_data.old$nps[[nps.intersect[[i]]]]$weight_ast[GO_data.old$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect], GO_data.new$nps[[nps.intersect[[i]]]]$weight_ast[GO_data.new$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect])
}

df <- data.frame(zeb = unlist(cor[[1]]), ast = unlist(cor[[2]]))


## Double for loop, for looping through nps/nts/synaptic/ion, then through each gene, then unlist to make ggplots

gene.lists <- c("nps", "nts", "synaptic", "ion")

cor <- list()
intersect <- list()
for (j in gene.lists) {
  cor[[j]] <- list()
  intersect[[j]] <- intersect(names(GO_data.old[[j]]), names(GO_data.new[[j]]))
  for (i in intersect[[j]]) {
    tf.intersect <- intersect(GO_data.old[[j]][[i]]$TF, GO_data.new[[j]][[i]]$TF)
    cor[[j]][[i]] <- data.frame(Zebrafish = cor(GO_data.old[[j]][[i]]$weight_zeb[GO_data.old[[j]][[i]]$TF %in% tf.intersect], GO_data.new[[j]][[i]]$weight_zeb[GO_data.new[[j]][[i]]$TF %in% tf.intersect]), Mexican_tetra = cor(GO_data.old[[j]][[i]]$weight_ast[GO_data.old[[j]][[i]]$TF %in% tf.intersect], GO_data.new[[j]][[i]]$weight_ast[GO_data.new[[j]][[i]]$TF %in% tf.intersect]))
  }
}
names(cor) <- c("Neuropeptide", "Neurotransmitter", "Synaptic", "Ion Channel")
# unlist each of the gene.lists
cor2 <- list()
for (i in 1:length(cor)) {
  cor2[[i]] <- Reduce(rbind, cor[[i]])
  cor2[[i]]$gene.name <- intersect[[i]]
  cor2[[i]] <- melt(cor2[[i]])
  cor2[[i]]$list <- names(cor)[[i]]
}

cor3 <- Reduce(rbind, cor2)
cor3$list <- factor(cor3$list, levels = c("Neuropeptide", "Neurotransmitter", "Synaptic", "Ion Channel"))

correlations.plot <- ggplot(cor3, aes(x = variable, y = value, colour = list)) + geom_boxplot(colour = "black", outlier.colour = "transparent") + geom_jitter(alpha = 0.5, size = 1) + facet_wrap(~list, scales = "free", nrow = 1) + theme_classic() + ylab("Pearson correlation") + xlab("") + theme(legend.position = "none", axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0))




## For each gene (np, nts, etc), I can plot the weight for zeb or ast in both old and new (many numbers, scatter plot)
## Can use this for comparing old and new vip or oxt plots
i <- grep("vip", nps.intersect)[[1]]
tf.intersect <- intersect(GO_data.old$nps[[nps.intersect[[i]]]]$TF, GO_data.new$nps[[nps.intersect[[i]]]]$TF)

df.gene <- data.frame(Zebrafish = GO_data.old$nps[[nps.intersect[[i]]]]$weight_zeb[GO_data.old$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect], Zebrafish_downsampled = GO_data.new$nps[[nps.intersect[[i]]]]$weight_zeb[GO_data.new$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect], Mexican_tetra = GO_data.old$nps[[nps.intersect[[i]]]]$weight_ast[GO_data.old$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect], Mexican_tetra_downsampled = GO_data.new$nps[[nps.intersect[[i]]]]$weight_ast[GO_data.new$nps[[nps.intersect[[i]]]]$TF %in% tf.intersect])
zeb.gene <- ggplot(df.gene, aes(x = Zebrafish, y = Zebrafish_downsampled)) + geom_point(size = 1) + theme_classic() + ylab("Downsampled weight") + xlab("Original weight") + ggtitle("Zebrafish")
ast.gene <- ggplot(df.gene, aes(x = Mexican_tetra, y = Mexican_tetra_downsampled)) + geom_point(size = 1) + theme_classic() + ylab("") + xlab("") + ggtitle("Mexican tetra")


# Make vip plots

vip.old <- plotPlotData(data = GO_data.old$nps, gene = "vip", quantile = 0.98, log = F) + theme(axis.text = element_text(size = 8), legend.position = "none") + ggtitle("vip") + ylab("RF weight - Zebrafish") + xlab("RF weight - Mexican tetra")
vip.new <- plotPlotData(data = GO_data.new$nps, gene = "vip", quantile = 0.98, log = F) + theme(axis.text = element_text(size = 8), axis.title = element_blank(), legend.position = "none") + ggtitle("vip", subtitle = "Downsampled")

# Put them all together

layout <- "
ABCD
EEEE"

dev.new()
vip.old + vip.new + zeb.gene + ast.gene + correlations.plot + plot_layout(design = layout, width = unit(c(25), "mm"), height = unit(c(25,40), "mm"))









