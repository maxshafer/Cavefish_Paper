library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(datapasta)
library(tidyr)
library(dendextend)
library(ggplot2)
library(cowplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration")
# load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")
# hypo <- hypo.integrated
# Subset to only Microglia

subsets <- readRDS("Hypo_integrated_130k_1500VFs_100Dims_subsets.rds")


## Make dendrogram of microglia cell types to see which are related

for (i in 1:length(subsets)) {
	Idents(subsets[[i]]) <- "integrated_SubclusterType"
}

subsets2 <- lapply(subsets, function(x) BuildClusterTree(x))
dend <- lapply(subsets2, function(x) as.dendrogram(Tool(x, slot = "BuildClusterTree")))
subset.dendrogram <- lapply(dend, function(x) ggplot(x %>% dendextend::set("labels_cex", value = 0.75), horiz = T, offset_labels = -2))


## Make tsne plots

cols <- c("lightgoldenrod1", "springgreen4","skyblue2")
cols2 <- c("#FF4500","#3E94D1")

subset.subcluster <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "integrated_SubclusterType", do.label = TRUE, no.legend = TRUE, do.return = TRUE, no.axes = T, pt.size = 2))

subset.species <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "species", do.label = FALSE, no.legend = T, do.return = TRUE, no.axes = T, pt.size = 2, cols = cols))

subset.sex <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "sex", do.label = FALSE, no.legend = T, do.return = TRUE, no.axes = T, pt.size = 2, cols = c("#FF4500", "#3E94D1")))

subset.species.2 <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "species.2", do.label = FALSE, no.legend = F, do.return = TRUE, no.axes = T, pt.size = 2, cols = cols))


# Make prop table for beside dendrogram

prop <- lapply(subsets, function(x) table(x@meta.data$integrated_SubclusterType, x@meta.data$species))
prop <- lapply(prop, function(x) as.data.frame(t(apply(x, 1, function(x) {x/sum(x)}))))

for (i in 1:length(prop)) {
	prop[[i]]$cell_type <- row.names(prop[[i]])
}

prop <- lapply(prop, function(x) x[,c(4,1,2,3)])
prop <- lapply(prop, function(x) x %>% gather(species_morph, freq, astyanax_cave:zebrafish))

for (i in 1:length(prop)) {
	prop[[i]]$cell_type <- factor(prop[[i]]$cell_type, levels = labels(dend[[i]]))
}

subset.prop <- list()
for (i in 1:length(prop)){
	subset.prop[[i]] <- ggplot(prop[[i]], aes(x=cell_type, y=freq, fill=species_morph)) + coord_flip() + geom_bar(stat="identity") + guides(fill=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank()) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols) + ylab("Species morph \n cluster frequency")
}

i <- 15
plot_grid(subset.dendrogram[[i]], subset.prop[[i]], subset.species[[i]], subset.subcluster[[i]], nrow = 1, rel_widths = c(1.5,.5,1,1))

FeaturePlot(subsets[[15]], features = c("txn"))
DotPlot(subsets[[15]], features = c("th", "th2", "prdx1", "txn"), split.by = "species", cols = c("blue", "red", "green"))

cols <- c("#FFBE00", "#FF4500","#3E94D1") #Red/Orange/Blue

# Load subcluster markers and separate!

markers.SubclusterType <- readRDS("Shafer_Hypo_markers.SubclusterType.rds")

markers.SubclusterType <- lapply(subtypes, function(x) markers.SubclusterType[grep(x, markers.SubclusterType$cluster),])

features.markers <- lapply(markers.SubclusterType, function(x) x %>% group_by(cluster) %>% top_n(5, avg_logFC))

for (i in 1:length(subsets)) {
	subset.features[[i]] <- DotPlot(object = subsets[[i]], genes.plot = unique(features.markers[[i]]$gene), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE, do.return = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}


plot_grids <- list()
for (i in 1:length(subtypes)) {
	plot_grids[[i]] <- plot_grid(subset.dendrogram[[i]], subset.prop[[i]], subset.subcluster[[i]], subset.features[[i]], ncol = 4, rel_widths = c(1,.2,.5,1))
}

png("Figures/AstMex_subcluster_analysis_figures.png", units = "in", res = 200, height = 100, width = 15)
plot_grid(plotlist = plot_grids, nrow = length(plot_grids))
dev.off()




# Average expression per cluster

normed.subcluster <- readRDS("AstMex_Hypo_normed.SubclusterType_matrix.rds")


scatter.plot <- 
ggplot(normed.subcluster, aes(x = GABA_1_4, y = GABA_1_3)) + geom_point(size = .5) 

+ xlab("log2(GABA_7_0 expression)") + ylab("log2(GABA_7_2 expression)")
+ geom_text_repel(aes(label = gene, colour = "darkred")) + guides(colour = FALSE)



data.to.plot <- log(cbind(norm.cluster[[1]], norm.cluster[[3]]))
names(data.to.plot) <- c("GABA_7_0", "GABA_7_2")
data.to.plot$gene <- rownames(data.to.plot)

genes.to.label <- c(subset.species.markers$gene[1:10], "hspa4a", "hspb1", "galn", "cort", "phactr2")

data.to.plot$gene <- ifelse(data.to.plot$gene %in% genes.to.label, data.to.plot$gene, "")

scatter.plot <- ggplot(data.to.plot, aes(x = GABA_7_0, y = GABA_7_2)) + geom_point(size = .5) + geom_text_repel(aes(label = gene, colour = "darkred")) + guides(colour = FALSE) + xlab("log2(GABA_7_0 expression)") + ylab("log2(GABA_7_2 expression)")


ggsave("Figures/AstMex_hypo_dot_plot_GABA_7.png", units = "in", dpi = 200, height = 6, width = 6, limitsize = FALSE)



# Plot ontology

clipr::write_clip(rownames(ciliated.markers.diff[ciliated.markers.diff$p_val_adj < 0.0001,]))
df_paste()

# Make vectors from pasted data

go.table <- data.frame(Term, Count, Percent, PValue, Benjamini)

go.plot <- 
ggplot(go.table, aes(x = Term, y = Percent, fill = Benjamini)) + geom_bar(stat="identity") + theme(axis.title.y = element_blank()) + coord_flip() + ggtitle("GO term biological process enrichment")





DimPlot(hypo, do.label = F, no.legend = T, no.axes = T, reduction.use = "FItSNE", do.return = F, vector.friendly = F, pt.size = 1, cells.highlight = (WhichCells(hypo, ident = "GABA_7_2")), cols.highlight = c("firebrick1"))


 
FeaturePlot(object = subset, no.axes = T, features.plot = c("hspb1"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = 2, reduction.use = "tsne")

DotPlot(object = subset, genes.plot = unique(gal.species.markers[gal.species.markers$p_val_adj < 1, "gene"]), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

DotPlot(object = subset, genes.plot = c("slc17a6a", "slc1a2b", "gad1b", "slc6a1b", "slc32a1", "rtn4rl2a", "rtn4rl2b", "meis2a", "cd9b", "hes6", "jdp2b", "rcan1a", "prelid3b", "jagn1a", "sat1b", "sat1a.2"), group.by = "orig.ident", plot.legend = TRUE, x.lab.rot = TRUE)

DotPlot(object = hypo, genes.plot = c("rtn4rl2a", "rtn4rl2b", "cd9b", "meis2a", "phactr3b", "cxcl14", "add3b", "tac1", "cart2"), group.by = "Subtype", plot.legend = TRUE, x.lab.rot = TRUE)





