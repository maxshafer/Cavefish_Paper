library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggdendro)
library(phylogram)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo")

# Load the dataset
load("AstMex_64k.Robj")
hypo <- hypo.ast
rm(hypo.ast)

cols0 <- c("lightgoldenrod1", "springgreen4")
cols1 <- c("firebrick1", "royalblue1") # For sex
cols2 <- c("lightgoldenrod3", "goldenrod3", "lightgoldenrod1", "goldenrod1", "lightgoldenrod2", "goldenrod2", "lightgoldenrod4", "goldenrod4", "springgreen1", "springgreen2", "springgreen3", "springgreen4", "palegreen1", "palegreen2", "palegreen3"," palegreen4") # Individal samples

cols3 <- c("goldenrod1", "goldenrod1", "lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1", "lightgoldenrod1", "darkorange1", "darkorange1", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4", "springgreen4"," springgreen4") # Individal caves/surface

cols4 <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788") # Good for orig ident, but might want similar colours

## Make TSNE graphs

png("Figures/hypo_cluster_plot_Subtype_nolab.png", height = 8, width = 8, units = "in", res = 500)
plot.Subtype <- DimPlot(object = hypo, group.by = "Subtype", reduction.use = "FItSNE", pt.size = .05, do.label = F,  label.size = 5, no.legend = TRUE, no.axes = T, vector.friendly = T)
dev.off()

png("Figures/hypo_cluster_plot_Subtype.png", height = 8, width = 8, units = "in", res = 500)
plot.Subtype <- DimPlot(object = hypo, group.by = "Subtype", reduction.use = "FItSNE", pt.size = .05, do.label = T,  label.size = 5, no.legend = TRUE, no.axes = T, vector.friendly = T)
dev.off()
png("Figures/hypo_cluster_plot_SubclusterType.png", height = 8, width = 8, units = "in", res = 500)
plot.SubclusterType <- DimPlot(object = hypo, group.by = "SubclusterType", reduction.use = "FItSNE", pt.size = .05, do.label = T,  label.size = 2, no.legend = TRUE, no.axes = T)
dev.off()
png("Figures/hypo_cluster_plot_orig_ident.png", height = 8, width = 8, units = "in", res = 500)
DimPlot(object = hypo, group.by = "orig.ident", reduction.use = "FItSNE", pt.size = .01, do.label = FALSE, no.legend = FALSE, no.axes = T, cols.use = cols2, vector.friendly = F) + theme(legend.position = c(0.91,0.87), legend.background = element_blank()) + guides(color = guide_legend(ncol = 2, override.aes = list(size = 5)))
dev.off()
png("Figures/hypo_cluster_plot_cave.png", height = 8, width = 8, units = "in", res = 500)
DimPlot(object = hypo, group.by = "orig.ident", reduction.use = "FItSNE", pt.size = .01, do.label = FALSE, no.legend = FALSE, no.axes = T, cols.use = cols3) + theme(legend.position = c(0.87,0.14), legend.background = element_blank()) + guides(color = guide_legend(ncol = 2, override.aes = list(size = 5)))
dev.off()
png("Figures/hypo_cluster_plot_species_morph.png", height = 8, width = 8, units = "in", res = 500)
DimPlot(object = hypo, group.by = "species", reduction.use = "FItSNE", do.label = FALSE, pt.size = .05, no.legend = FALSE, no.axes = T, cols.use = cols0) + theme(legend.position = c(0.9,0.9))
dev.off()
png("Figures/hypo_cluster_plot_sex.png", height = 8, width = 8, units = "in", res = 500) 
DimPlot(object = hypo, group.by = "sex", reduction.use = "FItSNE", do.label = FALSE, pt.size = .05, no.legend = FALSE, no.axes = T, cols.use = cols1) + theme(legend.position = c(0.9,0.9))
dev.off()

# Dot and Feature plots for major markers

png("Figures/hypo_featureplots_Subtype_markers.png", height = 4.5, width = 12, units = "in", res = 500) 
FeaturePlot(object = hypo, no.axes = T, features.plot = c("gng3", "slc17a6a", "gad1b", "her15.1", "prdx1", "otpb", "pfn1", "mpz", "mrc1a", "epd", "hopx", "hbaa2"), reduction.use = "FItSNE", cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, nCol = 6)
dev.off()

png("Figures/hypo_dotplots_Subtype_markers.png", height = 10, width = 6.5, units = "in", res = 500) 
DotPlot(object = hypo, genes.plot = rev(c("gng3", "slc17a6a", "gad1b", "her15.1", "prdx1", "otpb", "pfn1", "mpz", "mrc1a", "epd", "hopx", "hbaa2")), group.by = "Subtype", plot.legend = TRUE, x.lab.rot = TRUE) + theme(legend.position = "right")
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
prop.table[[3]] <- table(hypo@meta.data$Subtype, hypo@meta.data$species)
prop.table[[4]] <- table(hypo@meta.data$SubclusterType, hypo@meta.data$species)
prop.table[[5]] <- table(hypo@meta.data$Subtype, hypo@meta.data$sex)
prop.table[[6]] <- table(hypo@meta.data$SubclusterType, hypo@meta.data$sex)

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

prop.plots[[1]] <- prop.plots[[1]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols3)  + ylab("Sample Subtype frequency") + coord_flip() + guides(fill = guide_legend(title = "Origin (cave or surface)"))

prop.plots[[2]] <- prop.plots[[2]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols3) + ylab("Sample Subcluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Origin (cave or surface)"))

prop.plots[[3]] <- prop.plots[[3]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols0) + ylab("Species morph Subtype frequency") + coord_flip() + guides(fill = guide_legend(title = "Species morph"))

prop.plots[[4]] <- prop.plots[[4]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols0) + ylab("Species morph Subcluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Species morph"))

prop.plots[[5]] <- prop.plots[[5]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols) + ylab("Sex Subtype frequency") + coord_flip() + guides(fill = guide_legend(title = "Sex"))

prop.plots[[6]] <- prop.plots[[6]] + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols) + ylab("Sex Subcluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Sex"))

# Plot grid and save
png("Figures/AstMex_Hypo_Subtype_prop.png", units = "in", res = 250, height = 10, width = 20)
plot_grid(prop.plots[[1]], prop.plots[[3]], prop.plots[[5]], nrow = 1, rel_widths = c(1,1))
dev.off()

png("Figures/AstMex_Hypo_SubclusterType_prop.png", units = "in", res = 250, height = 25, width = 22.5)
plot_grid(prop.plots[[2]], prop.plots[[4]], prop.plots[[6]], nrow = 1, rel_widths = c(1,1))
dev.off()


# Get lists of peptides etc from zebrafish mine

go_lists <- list.files("/Volumes/Maxwell/R_Projects/AstMex_Hypo/")[grep("GO", list.files("/Volumes/Maxwell/R_Projects/AstMex_Hypo/"))]

go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))

go_lists <- lapply(go_lists, function(x) unique(x$V2))

test <- go_lists[[5]][go_lists[[5]] %in% row.names(hypo@raw.data)]
test <- markers.Subtype$gene[markers.Subtype$gene %in% go_lists[[4]]]

names(go_lists) <- list.files("/Volumes/Maxwell/R_Projects/AstMex_Hypo/")[grep("GO", list.files("/Volumes/Maxwell/R_Projects/AstMex_Hypo/"))]

names(go_lists)

# Plot
DotPlot(object = hypo, genes.plot = unique(test), group.by = "Subtype", plot.legend = TRUE, x.lab.rot = TRUE) + theme(legend.position = "right")


