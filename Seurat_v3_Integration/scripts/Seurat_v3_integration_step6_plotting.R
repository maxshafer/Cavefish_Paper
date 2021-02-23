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

png("Figures/Hypo_integrated_umap_Cluster.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "umap", group.by = "integrated_Cluster", pt.size = .1, label = TRUE, shuffle = T) + NoAxes() + NoLegend()
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

marker.sub <- DotPlot(hypo.integrated, features = rev(unique(unlist(genes.to.plot))), group.by = "integrated_Cluster", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1), axis.title = element_blank())
marker.sub <- marker.sub + theme(axis.text = element_text(size = 6))

png("Figures/Hypo_integrated_dotplots_Cluster_markers.png", height = 10, width = 6.5, units = "in", res = 500) 
marker.sub
dev.off()

# Make neurotransmitter dot plot for Neuronal subclusters
Idents(hypo.integrated) <- "integrated_Cluster"
hypo.integrated.neuronal <- subset(hypo.integrated, ident = levels(Idents(hypo.integrated))[grep("Neuronal", levels(Idents(hypo.integrated)))])

nt.plot <- DotPlot(object = hypo.integrated.neuronal, features = rev(c("slc17a6a", "slc17a6b", "gad1b", "gad2")), group.by = "integrated_Subcluster") + theme(legend.position = "right") + RotatedAxis() + scale_color_viridis()
nt.plot <- nt.plot + theme(axis.text = element_text(size = 6), axis.title = element_blank())

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
props$Var1 <- factor(props$Var1, levels = names(sort(type_total)))

prop.plot <- ggplot(props, aes(x = value, y = Var1, width = (type_total), fill = Var2)) + geom_bar(stat = "identity", position = "fill", colour = "black", size = 0) + scale_fill_manual(values = cols) + theme(panel.spacing.y = unit(.25, "mm"), strip.text.y = element_blank(), strip.background = element_blank(), axis.ticks.length.y =  unit(.5, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10), legend.title = element_blank()) + xlab("Prop. of total species cells") + ylab(element_blank()) + geom_vline(xintercept = c(0.3333, 0.6666), linetype = "dashed", color = "black") + facet_grid(rows = vars(Var1), scales = "free", space = "free")
prop.plot <- prop.plot + theme(axis.title.x = element_blank())
pdf("Figures/Hypo_integrated_Cluster_barplots.pdf", width = 10, height = 3)
prop.plot
dev.off()

## Plot plots using patchwork
# Figure S4 - marker genes and proportion

marker.sub + prop.plot + plot_layout(height = unit(c(150), c("mm")), width = unit(c(75,20), c("mm")), guides = "collect")

# Figure S5 - nt usage

nt.plot + plot_layour()
