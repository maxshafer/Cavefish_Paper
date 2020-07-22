library(Seurat)
library(Matrix)
library(dplyr)
library(datapasta)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/")

load("Shafer_Hypo_66k.Robj")

# Subset to only Microglia

hypo <- SetAllIdent(hypo, id = "Subtype")

subtypes <- levels(hypo@ident)

set.seed(1188)

processSubset <- function(object = object, ident = ident, perplexity = 30) {
	object <- SetAllIdent(object, id = "Subtype")
	subset <- SubsetData(hypo, ident.use = ident)
	subset <- FindVariableGenes(subset, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = FALSE)
	subset <- ScaleData(subset, genes.use = subset@var.genes, do.par = TRUE, num.cores = 6)
	subset <- RunPCA(subset, pc.genes = subset@var.genes, pcs.compute = 15)
	subset <- RunTSNE(subset, reduction.use = "pca", dims.use = 1:15, nthreads = 6, do.fast = TRUE, check_duplicates = FALSE, perplexity = perplexity)
	return(subset)
}

subsets <- lapply(subtypes, function(x) processSubset(object = hypo, id = x))
names(subsets) <- subtypes
subsets[[10]] <- processSubset(hypo, ident = subtypes[10], perplexity = 30)

## Make dendrogram of microglia cell types to see which are related

subsets <- lapply(subsets, function(x) SetAllIdent(x, id = "SubclusterType"))
subsets[c(1:10,12,15,19:27,30:31,33:36,40:41)] <- lapply(subsets[c(1:10,12,15,19:27,30:31,33:36,40:41)], function(x) BuildClusterTree(x))
dend <- lapply(subsets, function(x) as.dendrogram(x@cluster.tree[[1]]))
subset.dendrogram <- lapply(dend, function(x) ggplot(x %>% dendextend::set("labels_cex", value = 0.75), horiz = T, offset_labels = -2))

subset.dendrogram <- lapply(subset.dendrogram, function(x) x + theme(plot.margin = unit(c(0,0,0,0), units = "cm")))


# Save subsets for re-loading

saveRDS(subsets, file = "Shafer_Hypo_subtype_subsets.rds")

subsets <- readRDS("AstMex_Hypo_subtype_subsets.rds")

## Make tsne plots

subset.subcluster <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "SubclusterType", do.label = TRUE, no.legend = TRUE, do.return = TRUE, no.axes = T, pt.size = 2))

subset.sex <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "sex", do.label = FALSE, no.legend = T, do.return = TRUE, no.axes = T, pt.size = 2, cols.use = c("indianred1", "firebrick1", "royalblue1", "skyblue1")))


# Make prop table for beside dendrogram

prop <- lapply(subsets, function(x) table(x@meta.data$SubclusterType, x@meta.data$sex))
prop <- lapply(prop, function(x) as.data.frame(t(apply(x, 1, function(x) {x/sum(x)}))))

for (i in 1:length(prop)) {
	prop[[i]]$cell_type <- row.names(prop[[i]])
}

prop <- lapply(prop, function(x) x[,c(5,1,2,3,4)])
prop <- lapply(prop, function(x) x %>% gather(species_morph, freq, female1:male2))

for (i in 1:length(prop)) {
	prop[[i]]$cell_type <- factor(prop[[i]]$cell_type, levels = labels(dend[[i]]))
}

subset.prop <- list()
for (i in 1:length(prop)){
	subset.prop[[i]] <- ggplot(prop[[i]], aes(x=cell_type, y=freq, fill=species_morph)) + coord_flip() + geom_bar(stat="identity") + guides(fill=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank()) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = c("indianred1", "firebrick1", "royalblue1", "skyblue1")) + ylab("Species morph \n cluster frequency")
}





cols <- c("#FFBE00", "#FF4500","#3E94D1") #Red/Orange/Blue

# Load subcluster markers and separate!

markers.SubclusterType <- readRDS("Shafer_Hypo_markers.SubclusterType.rds")

features.markers <- markers.SubclusterType %>% group_by(SubclusterType) %>% top_n(-5, value)

subset.features <- list()
for (i in 1:length(subsets)) {
	if (length(features.markers[features.markers$L1 == subtypes[[i]], "gene"]$gene) > 0){
			subset.features[[i]] <- DotPlot(object = subsets[[i]], genes.plot = unique(features.markers[features.markers$L1 == subtypes[[i]], "gene"]$gene), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE, do.return = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

	}
}


## Plot the grids!

plot_grids <- list()
for (i in 1:length(subtypes)) {
	plot_grids[[i]] <- plot_grid(subset.dendrogram[[i]], subset.prop[[i]], subset.subcluster[[i]], subset.features[[i]], ncol = 4, rel_widths = c(1,.2,.5,1))
}

png("Figures/AstMex_subcluster_analysis_figures.png", units = "in", res = 200, height = 100, width = 15)
plot_grid(plotlist = plot_grids, nrow = length(plot_grids))
dev.off()


# Average expression per cluster

normCluster <- function(object = hypo, idents = "Subtype") {
	object <- SetAllIdent(object, id = as.character(idents))
	subsets <- lapply(levels(object@ident), function(x) SubsetData(object, ident.use = x))
	norm.cluster <- lapply(subsets, function(x) data.frame(mean.exp = apply(expm1(x@data), 1, function(x) mean(x))))
	names(norm.cluster) <- levels(object@ident)
	return(norm.cluster)iden
}

id.types <- c("Subtype", "SubclusterType", "species_subtype", "species_subcluster")

norm.cluster.data <- lapply(id.types[1:2], function(x) normCluster(object = hypo, idents = x))

cluster1 <- "GABA_4_0"
cluster2 <- "GABA_4_1"

hypo <- SetAllIdent(hypo, id = "SubclusterType")
test <- FindMarkers(hypo, ident.1 = cluster1, ident.2 = cluster2)

genes.to.label <- c(norm.cluster.data[[2]]$gene[1:10], c(row.names(test[1:20,])))

data.to.plot <- log(cbind(norm.cluster.data[[2]][[cluster1]], norm.cluster.data[[2]][[cluster2]]))
names(data.to.plot) <- c(cluster1, cluster2)
data.to.plot$gene <- rownames(data.to.plot)


data.to.plot$gene <- ifelse(data.to.plot$gene %in% genes.to.label, data.to.plot$gene, "")

scatter.plot <- ggplot(data.to.plot, aes(x = GABA_4_0, y = GABA_4_1)) + geom_point(size = .5) + geom_text_repel(aes(label = gene, colour = "darkred")) + guides(colour = FALSE) + xlab(paste("log2(", cluster1, " expression)", sep = "")) + ylab(paste("log2(", cluster2, " expression)", sep = ""))

scatter.plot.GABA_4 <- scatter.plot
scatter.plot.Glut_6 <- scatter.plot

i <- 25
plot_grid(plot_grid(plot_grid(subset.dendrogram[[i]], subset.prop[[i]], subset.sex[[i]], subset.subcluster[[i]], ncol = 4, rel_widths = c(1,.2,.5,.5), labels = c("A", "", "B", "C")), plot_grid(subset.features[[i]], scatter.plot.Glut_6, ncol = 2, rel_widths = c(2.5,1), labels = c("D", "E")), nrow = 2))

ggsave("Figures/Shafer_hypo_cluster_plots_GABA_4.png", units = "in", dpi = 200, height = 10, width = 20, limitsize = FALSE)

i <- 8
plot_grid(plot_grid(plot_grid(subset.dendrogram[[i]], subset.prop[[i]], subset.sex[[i]], subset.subcluster[[i]], ncol = 4, rel_widths = c(1,.2,.5,.5), labels = c("A", "", "B", "C")), plot_grid(subset.features[[i]], scatter.plot.GABA_4, ncol = 2, rel_widths = c(2.5,1), labels = c("D", "E")), nrow = 2))

ggsave("Figures/Shafer_hypo_cluster_plots_GABA_4.png", units = "in", dpi = 200, height = 10, width = 20, limitsize = FALSE)












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





