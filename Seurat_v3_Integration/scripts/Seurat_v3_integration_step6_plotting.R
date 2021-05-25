library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(viridis)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

# Save plots (this doesn't work in for loop, dunno why not)
cols <- c("#FDE725FF", "#22A884FF", "#414487FF")

cols2 <- c("#FDE725FF", "#414487FF")
cols3 <- c("#FDE725FF", "#22A884FF", "#414487FF")

# Save Species, species.2, Cluster and subcluster plots + UMAPs
png("Figures/Hypo_integrated_tsne_batch_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "basetsne", group.by = "species", pt.size = .1, label = F, cols = cols, shuffle = T) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_umap_batch_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "baseumap", group.by = "species", pt.size = .1, label = F, cols = cols, shuffle = T) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "species", pt.size = .1, label = F, cols = cols, shuffle = T) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_species2.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "species.2", pt.size = .1, label = F, cols = cols2, shuffle = T) + NoAxes() + theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_Cluster.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Cluster", pt.size = .1, label = TRUE, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_umap_species2.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "umap", group.by = "species.2", pt.size = .1, label = F, cols = cols2, shuffle = T) + NoAxes() + theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_umap_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "umap", group.by = "species", pt.size = .1, label = F, cols = cols, shuffle = T) + NoAxes() + theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_umap_Cluster.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "umap", group.by = "integrated_Cluster", pt.size = .1, label = TRUE, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_umap_Cluster_nolabel.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "umap", group.by = "integrated_Cluster", pt.size = .1, label = F, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_Cluster_nolabel.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Cluster", pt.size = .1, label = F, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_Subcluster_nolabel.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subcluster", pt.size = .1, label = F, shuffle = T) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_Subcluster.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subcluster", pt.size = .1, label = T) + NoAxes() + NoLegend()
dev.off()


# DotPlots for major markers

png("Figures/Hypo_integrated_dotplots_markers.png", height = 10, width = 6.5, units = "in", res = 500) 
DotPlot(object = hypo.integrated, features = rev(c("pfn1", "mrc1a","gng3", "slc17a6a", "gad1b", "prdx1","mpz", "her15.1", "epd","arl13b", "hbaa2", "hopx")), group.by = "integrated_Cluster") + theme(legend.position = "right") + RotatedAxis() + scale_color_viridis()
dev.off()

# DotPlots for Cluster marker genes

gene.lists <- readRDS("drift_gene_lists_2.rds")

genes.to.plot <- lapply(gene.lists[[1]], function(x) row.names(x)[1:5])

# Replace markers for Leucocytes with other genes (#1 is good)
genes.to.plot$Leucocytes <- row.names(gene.lists[[1]]$Leucocytes)[c(1,6,9,16,57,60)]

marker.sub <- DotPlot(hypo.integrated, features = rev(unique(unlist(genes.to.plot))), group.by = "integrated_Cluster", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1), axis.title = element_blank())

### Make a prop plot with x-axis scaled by cluster size

props <- table(hypo.integrated@meta.data$integrated_Cluster, hypo.integrated@meta.data$species)
type_total <- apply(props, 1, function(x) sum(x[1],x[2],x[3]))
type_total <- (type_total/sum(type_total)*1000)
props[,1] <- props[,1]/sum(props[,1])*100
props[,2] <- props[,2]/sum(props[,2])*100
props[,3] <- props[,3]/sum(props[,3])*100
props <- melt(props)
props$type_total <- rep(type_total, 3)

# Order by total population level (can change to another thing)
# props$Var1 <- factor(props$Var1, levels = names(sort(type_total)))

prop.plot <- ggplot(props, aes(x = value, y = Var1, width = (type_total), fill = Var2)) + geom_bar(stat = "identity", position = "fill", colour = "black", size = 0) + scale_fill_manual(values = cols) + theme(panel.spacing.y = unit(.25, "mm"), strip.text.y = element_blank(), strip.background = element_blank(), axis.ticks.length.y =  unit(.5, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10), legend.title = element_blank()) + xlab("Prop. of total species cells") + ylab(element_blank()) + geom_vline(xintercept = c(0.3333, 0.6666), linetype = "dashed", color = "black") + facet_grid(rows = vars(Var1), scales = "free", space = "free")
prop.plot <- prop.plot + theme(axis.title.x = element_blank())

# Make % species cell histogram

props <- table(hypo.integrated@meta.data$integrated_Subcluster, hypo.integrated@meta.data$species.2)
colnames(props) <- c("Mexican tetra", "Zebrafish")
props <- props/rowSums(props)
props <- melt(props)
props$Var2 <- factor(props$Var2, levels = c("Zebrafish", "Mexican tetra")) 

density.plot <- ggplot(props[props$value > 0.5,], aes(x = value, group = Var2, fill = Var2, colour = Var2)) + geom_density(alpha = 0.25) + theme_classic() + scale_colour_viridis_d() + scale_fill_viridis_d() + theme(legend.position = c(0.65,0.75), legend.title = element_blank())
density.plot <- density.plot + ylab("Density") + xlab("Fraction of subcluster") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10))
## Plot plots using patchwork
# Figure S4 - marker genes and proportion

pdf("Figures/Hypo_integrated_Figure-S3.pdf", height = 12, width = 8) 
marker.sub + prop.plot + density.plot + plot_layout(ncol = 2, height = unit(c(200, 20), c("mm")), width = unit(c(80,22,80), c("mm")), guides = "collect")
dev.off()

