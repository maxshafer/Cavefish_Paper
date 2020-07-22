library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/")

# Load object
load("Shafer_Hypo_66k.Robj")
hypo <- hypo.zeb
rm(hypo.zeb)

# cols <- c("lightgoldenrod1", "springgreen4")
cols <- c("indianred1", "firebrick1", "royalblue1", "skyblue1")
cols2 <- c("royalblue1", "skyblue2", "skyblue3", "skyblue4", "firebrick1", "firebrick2", "firebrick3", "firebrick4", "royalblue2", "royalblue3", "royalblue4", "indianred1", "indianred2", "indianred3", "indianred"," skyblue1")
cols3 <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")

## Make TSNE graphs

## For figures

png("Figures/hypo_cluster_plot_Subtype_nolab.png", height = 8, width = 8, units = "in", res = 500)
plot.Subtype <- DimPlot(object = hypo, group.by = "Subtype", reduction.use = "tsne", pt.size = .05, do.label = F,  label.size = 5, no.legend = TRUE, no.axes = T, vector.friendly = T)
dev.off()

pdf("Figures/hypo_cluster_plot_Subtype.pdf", height = 8, width = 8)
plot.Subtype <- DimPlot(object = hypo, group.by = "Subtype", reduction.use = "tsne", pt.size = .05, do.label = T,  label.size = 5, no.legend = TRUE, no.axes = T, vector.friendly = T)
dev.off()
png("Figures/hypo_cluster_plot_SubclusterType.png", height = 8, width = 8, units = "in", res = 500)
plot.SubclusterType <- TSNEPlot(object = hypo, group.by = "SubclusterType", pt.size = .05, do.label = T, label.size = 2, no.legend = TRUE, no.axes = T)
dev.off()
png("Figures/hypo_cluster_plot_orig_ident.png", height = 8, width = 8, units = "in", res = 500)
TSNEPlot(object = hypo, group.by = "orig.ident", pt.size = .01, do.label = FALSE, no.legend = FALSE, no.axes = T, colors.use = cols3) + theme(legend.position = c(0.9,0.85), legend.background = element_blank()) + guides(color = guide_legend(ncol = 2, override.aes = list(size = 5)))
dev.off()
png("Figures/hypo_cluster_plot_sex.png", height = 8, width = 8, units = "in", res = 500) 
TSNEPlot(object = hypo, group.by = "sex", do.label = FALSE, pt.size = .05, no.legend = FALSE, no.axes = T, colors.use = cols) + theme(legend.position = c(0.9,0.9))
dev.off()

png("Figures/hypo_featureplots_Subtype_markers.png", height = 4.5, width = 12, units = "in", res = 500) 
FeaturePlot(object = hypo, no.axes = T, features.plot = c("gng3", "slc17a6b", "gad2", "her4.2", "prdx1", "otpa", "pfn1", "olig2", "mrc1a", "epd", "hopx", "ba1"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, nCol = 6)
dev.off()

png("Figures/hypo_dotplots_Subtype_markers.png", height = 10, width = 6.5, units = "in", res = 500) 
DotPlot(object = hypo, genes.plot = rev(c("gng3", "slc17a6b", "gad2", "her4.2", "prdx1", "otpa", "pfn1", "olig2", "mrc1a", "epd", "hopx", "ba1")), plot.legend = TRUE, x.lab.rot = TRUE) + theme(legend.position = "right")
dev.off()


# Proportion figures for sex and orig.ident

ids <- c("Subtype", "SubclusterType")

my_colour_palette <- list()
for (i in 1:length(ids)) {
	hypo <- SetAllIdent(hypo, id = ids[[i]])
	colours <- hue_pal()(length(levels(hypo@ident)))
	names(colours) <- levels(hypo@ident)
	my_colour_palette[[i]] <- colours
}
my_colour_palette[[3]] <- my_colour_palette[[1]]
my_colour_palette[[4]] <- my_colour_palette[[2]]

# Make tables of cell type proportions
prop.table <- list()
prop.table[[1]] <- table(hypo@meta.data$Subtype, hypo@meta.data$orig.ident)
prop.table[[2]] <- table(hypo@meta.data$SubclusterType, hypo@meta.data$orig.ident)
prop.table[[3]] <- table(hypo@meta.data$Subtype, hypo@meta.data$sex)
prop.table[[4]] <- table(hypo@meta.data$SubclusterType, hypo@meta.data$sex)

prop.table <- lapply(prop.table, function(x) as.data.frame(t(apply(x, 1, function(y) {y/sum(y)}))))

for (i in 1:length(prop.table)) {
	prop.table[[i]]$cell_type <- row.names(prop.table[[i]])
	# prop.table[[i]] <- prop.table[[i]][,c(5,1,2,3,4)]
	prop.table[[i]] <- melt(prop.table[[i]])
	prop.table[[i]]$cell_type <- as.factor(prop.table[[i]]$cell_type)
}

prop.plots <- list()
for (i in 1:length(prop.table)) {
	prop.plots[[i]] <- ggplot(prop.table[[i]], aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity")
}

prop.plots[[1]] <- prop.plots[[1]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols3)  + ylab("Sample Subtype frequency") + coord_flip() + guides(fill = guide_legend(title = "10x Sample"))

prop.plots[[2]] <- prop.plots[[2]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols3) + ylab("Sample Subcluster frequency") + coord_flip() + guides(fill = guide_legend(title = "10x Sample"))


prop.plots[[3]] <- prop.plots[[3]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols) + ylab("Sex Subtype frequency") + coord_flip() + guides(fill = guide_legend(title = "Sex sample"))

prop.plots[[4]] <- prop.plots[[4]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols) + ylab("Sex Subcluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Sex sample"))


# Plot grid and save
png("Figures/Shafer_Hypo_Subtype_prop.png", units = "in", res = 250, height = 10, width = 13)
plot_grid(prop.plots[[1]], prop.plots[[3]], nrow = 1, rel_widths = c(1,1))
dev.off()

png("Figures/Shafer_Hypo_SubclusterType_prop.png", units = "in", res = 250, height = 25, width = 15)
plot_grid(prop.plots[[2]], prop.plots[[4]], nrow = 1, rel_widths = c(1,1))
dev.off()



# Subcluster Figure graphs

png("Figures/hypo_dotplots_hcrt.png", units = "in", res = 250, height = 20, width = 6)
DotPlot(object = hypo, genes.plot = rev(c("lhx9", "npvf", "NP5", "hcrt")), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE) + theme(legend.position = "right")
dev.off()

png("Figures/hypo_cluster_plot_highlight_Glut_18_3.png", height = 8, width = 8, units = "in", res = 250)
DimPlot(hypo, do.label = F, no.legend = T, no.axes = T, reduction.use = "tsne", do.return = F, vector.friendly = F, pt.size = 1, cells.highlight = WhichCells(hypo, ident = c("Glut_6_3")), cols.highlight = c("firebrick1"))
dev.off()

hypo <- SetAllIdent(hypo, id = "Subtype")
subset <- SubsetData(hypo, ident.use = "Glut_6")

subset <- NormalizeData(object = subset, normalization.method = "LogNormalize", scale.factor = 10000)
subset <- FindVariableGenes(object = subset, mean.function = ExpMean, dispersion.function = LogVMR, do.plot = TRUE)
subset <- ScaleData(object = subset, genes.use = subset@var.genes, do.par = TRUE, num.cores = 6)
length(x = subset@var.genes)

subset <- RunPCA(object = subset, pc.genes = subset@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)
subset <- RunTSNE(object = subset, reduction.use = "pca", dims.use = 1:5, nthreads = 6, do.fast = TRUE, perplexity = 10)

DimPlot(subset, do.label = T, no.legend = T, no.axes = T, group.by = "SubclusterType", reduction.use = "tsne", do.return = F, vector.friendly = F, pt.size = 2)

FeaturePlot(object = hypo, no.axes = T, features.plot = c("rfx4"), cols.use = c("grey85", "blue"), pt.size = 2, no.legend = T, min.cutoff = "q9")

DotPlot(object = subset, genes.plot = unique(markers.SubclusterType$gene[markers.SubclusterType$L1 == "Glut_6"]), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE) + theme(legend.position = "right")


FeaturePlot(object = subset, no.axes = T, features.plot = c("npvf", "hcrt"), overlay = T, cols.use = c("grey85", "blue", "red", "green"), pt.size = 2, no.legend = F)
dev.new()
FeaturePlot(object = subset, no.axes = T, features.plot = c("npvf", "NP5"), overlay = T, cols.use = c("grey85", "blue", "red", "green"), pt.size = 2, no.legend = F)
dev.new()
FeaturePlot(object = subset, no.axes = T, features.plot = c("hcrt", "NP5"), overlay = T, cols.use = c("grey85", "blue", "red", "green"), pt.size = 2, no.legend = F)


# Get lists of peptides etc from zebrafish mine

go_lists <- list.files()[grep("GO", list.files())]

go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))

go_lists <- lapply(go_lists, function(x) unique(x$V2))
go_lists <- lapply(go_lists, function(x) x[x %in% markers.SubclusterType$gene])

names(go_lists) <- list.files()[grep("GO", list.files())]

names(go_lists)

# Plot
DotPlot(object = hypo, genes.plot = unique(go_lists[[5]]), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE, do.return = T) + theme(legend.position = "right", axis.text.x = element_text(hjust = 1)) + coord_flip()

ggsave("Figures/Shafer_hypo_SubclusterType_peptides.png", units = "in", width = 25, height = 10)




