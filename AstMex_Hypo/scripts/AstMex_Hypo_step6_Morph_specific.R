library(Seurat)
library(Matrix)
library(dplyr)
library(datapasta)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/")

load("AstMex_64k.Robj")

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

subsets <- lapply(subtypes[31:36], function(x) processSubset(object = hypo, id = x))

subsets[[10]] <- processSubset(hypo, ident = subtypes[10], perplexity = 30)

## Make dendrogram of microglia cell types to see which are related

subsets <- lapply(subsets, function(x) SetAllIdent(x, id = "SubclusterType"))
subsets[c(1:9,11:22,25:29,31:36)] <- lapply(subsets[c(1:9,11:22,25:29,31:36)], function(x) BuildClusterTree(x))
dend <- lapply(subsets, function(x) as.dendrogram(x@cluster.tree[[1]]))
subset.dendrogram <- lapply(dend, function(x) ggplot(x %>% dendextend::set("labels_cex", value = 0.75), horiz = T, offset_labels = -2))


# Save subsets for re-loading

saveRDS(subsets, file = "AstMex_Hypo_subtype_subsets.rds")

subsets <- readRDS("AstMex_Hypo_subtype_subsets.rds")

## Make tsne plots

subset.subcluster <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "SubclusterType", do.label = TRUE, no.legend = TRUE, do.return = TRUE, no.axes = T, pt.size = 2))

subset.species <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "species", do.label = FALSE, no.legend = T, do.return = TRUE, no.axes = T, pt.size = 2, cols.use = c("lightgoldenrod2", "springgreen4")))

subset.sex <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "sex", do.label = FALSE, no.legend = T, do.return = TRUE, no.axes = T, pt.size = 2, cols.use = c("#FF4500", "#3E94D1")))

subset.morph <- lapply(subsets, function(x) DimPlot(x, reduction.use = "tsne", group.by = "orig.ident", do.label = FALSE, no.legend = F, do.return = TRUE, no.axes = T, pt.size = 2, cols.use = c("gold", "gold1","firebrick", "firebrick1", "firebrick2", "firebrick3","orange", "orange3", "skyblue", "skyblue1", "skyblue2", "skyblue3", "aquamarine", "aquamarine1", "aquamarine2", "aquamarine3")))


# Make prop table for beside dendrogram

prop <- lapply(subsets, function(x) table(x@meta.data$SubclusterType, x@meta.data$species))
prop <- lapply(prop, function(x) as.data.frame(t(apply(x, 1, function(x) {x/sum(x)}))))

for (i in 1:length(prop)) {
	prop[[i]]$cell_type <- row.names(prop[[i]])
}

prop <- lapply(prop, function(x) x[,c(3,1,2)])
prop <- lapply(prop, function(x) x %>% gather(species_morph, freq, astyanax_cave:astyanax_surface))

for (i in 1:length(prop)) {
	prop[[i]]$cell_type <- factor(prop[[i]]$cell_type, levels = labels(dend[[i]]))
}

subset.prop <- list()
for (i in 1:length(prop)){
	subset.prop[[i]] <- ggplot(prop[[i]], aes(x=cell_type, y=freq, fill=species_morph)) + coord_flip() + geom_bar(stat="identity") + guides(fill=FALSE) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank()) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = c("lightgoldenrod1", "springgreen4")) + ylab("Species morph \n cluster frequency")
}





cols <- c("#FFBE00", "#FF4500","#3E94D1") #Red/Orange/Blue

# Load subcluster markers and separate!

markers.SubclusterType <- readRDS("Shafer_Hypo_markers.SubclusterType.rds")

markers.SubclusterType <- lapply(subtypes, function(x) markers.SubclusterType[grep(x, markers.SubclusterType$cluster),])

features.markers <- lapply(markers.SubclusterType, function(x) x %>% group_by(cluster) %>% top_n(5, avg_logFC))

subset.features <- list()
for (i in 1:length(subsets)) {
	subset.features[[i]] <- DotPlot(object = subsets[[i]], genes.plot = unique(features.markers[[i]]$gene), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE, do.return = TRUE) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
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

norm.cluster.data[[4]] <- lapply(id.types[[4]], function(x) normCluster(object = hypo, idents = x))

cluster1 <- "GABA_4_0"
cluster2 <- "GABA_4_1"

hypo <- SetAllIdent(hypo, id = "SubclusterType")
test <- FindMarkers(hypo, ident.1 = cluster1, ident.2 = cluster2)

genes.to.label <- c(norm.cluster.data[[2]]$gene[1:10], c(row.names(test[1:50,])))

data.to.plot <- log(cbind(norm.cluster.subclustertype[[cluster1]], norm.cluster.subclustertype[[cluster2]]))
names(data.to.plot) <- c(cluster1, cluster2)
data.to.plot$gene <- rownames(data.to.plot)


data.to.plot$gene <- ifelse(data.to.plot$gene %in% genes.to.label, data.to.plot$gene, "")

scatter.plot <- ggplot(data.to.plot, aes(x = GABA_4_0, y = GABA_4_1)) + geom_point(size = .5) + geom_text_repel(aes(label = gene, colour = "darkred")) + guides(colour = FALSE) + xlab(paste("log2(", cluster1, " expression)", sep = "")) + ylab(paste("log2(", cluster2, " expression)", sep = ""))


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





