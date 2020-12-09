library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(viridis)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v3.rds")

# Save plots (this doesn't work in for loop, dunno why not)
cols <- c("#FDE725FF", "#22A884FF", "#414487FF")

cols2 <- c("#FDE725FF", "#414487FF")
cols3 <- c("#FDE725FF", "#22A884FF", "#414487FF")

# Save Species, species.2, subtype and subcluster plots
png("Figures/Hypo_integrated_tsne_batch_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "basetsne", group.by = "species", pt.size = .1, label = F, cols = cols, shuffle = T) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "species", pt.size = .1, label = F, cols = cols, shuffle = T) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_species2.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "species.2", pt.size = .1, label = F, cols = cols2, shuffle = T) + NoAxes() + theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_subtype.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subtype", pt.size = .1, label = TRUE, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_subtype_nolabel.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subtype", pt.size = .1, label = F, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_subclustertype_nolabel.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_SubclusterType", pt.size = .1, label = F, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_subclustertype.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_SubclusterType", pt.size = .1, label = T) + NoAxes() + NoLegend()
dev.off()


# DotPlots for major markers

png("Figures/Hypo_integrated_dotplots_markers.png", height = 10, width = 6.5, units = "in", res = 500) 
DotPlot(object = hypo.integrated, features = rev(c("gng3", "slc17a6a", "gad1b", "her15.1", "prdx1", "otpb", "pfn1", "mpz", "mrc1a", "epd", "hopx", "hbaa2")), group.by = "Subtype") + theme(legend.position = "right") + RotatedAxis() + scale_color_viridis()
dev.off()

# DotPlots for Subtype marker genes

gene.lists <- readRDS("drift_gene_lists.rds")

genes.to.plot <- lapply(gene.lists[[1]], function(x) row.names(x)[1:5])

marker.sub <- DotPlot(hypo.integrated, features = rev(unique(unlist(genes.to.plot))), group.by = "integrated_Subtype", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1), axis.title = element_blank())

png("Figures/Hypo_integrated_dotplots_Subtype_markers.png", height = 10, width = 6.5, units = "in", res = 500) 
marker.sub
dev.off()


### Make a prop plot with x-axis scaled by cluster size

props <- table(hypo.integrated@meta.data$integrated_Subtype, hypo.integrated@meta.data$species)
type_total <- apply(props, 1, function(x) sum(x[1],x[2],x[3]))
type_total <- (type_total/sum(type_total)*1000)
props[,1] <- props[,1]/sum(props[,1])*100
props[,2] <- props[,2]/sum(props[,2])*100
props[,3] <- props[,3]/sum(props[,3])*100
props <- melt(props)
props$type_total <- rep(type_total, 3)

# Order by total population level (can change to another thing)
props$Var1 <- factor(props$Var1, levels = names(sort(type_total)))

prop.plot <- ggplot(props, aes(x = value, y = Var1, width = (type_total), fill = Var2)) + geom_bar(stat = "identity", position = "fill", colour = "black", size = 0) + scale_fill_manual(values = cols) + theme(panel.spacing.y = unit(.25, "mm"), strip.text.y = element_blank(), strip.background = element_blank(), axis.ticks.length.y =  unit(.5, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10), legend.title = element_blank()) + xlab("Prop. of total species cells") + ylab(element_blank()) + geom_vline(xintercept = c(0.3333, 0.6666), linetype = "dashed", color = "black") + facet_grid(rows = vars(Var1), scales = "free", space = "free")

pdf("Figures/Hypo_integrated_Subtype_barplots.pdf", width = 10, height = 3)
prop.plot
dev.off()

library(patchwork)

marker.sub + prop.plot + plot_layout(widths = c(3,1), guides = "collect")

library(scales)

highlight.cells <- lapply(levels(hypo.integrated@meta.data$integrated_Subtype)[grep("GABA", levels(hypo.integrated@meta.data$integrated_Subtype))], function(x) WhichCells(hypo.integrated, idents = x))
colors <- hue_pal()(length(highlight.cells))

DimPlot(hypo.integrated, reduction = "tsne", cells.highlight = highlight.cells, cols.highlight = colors, sizes.highlight = 0.1)
