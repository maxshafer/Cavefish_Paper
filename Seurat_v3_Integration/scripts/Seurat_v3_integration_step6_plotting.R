library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(reshape2)
library(ggdendro)
library(phylogram)
library(ggtree)
library(grid)
library(ggmap)
library(dendextend)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration")
load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

subsets.all <- readRDS("Hypo_integrated_130k_1500VFs_100Dims_subsets.rds")

# Save plots (this doesn't work in for loop, dunno why not)
cols <- c("#FDE725FF", "#22A884FF", "#414487FF")

cols2 <- c("#FDE725FF", "#414487FF")
cols3 <- c("#FDE725FF", "#22A884FF", "#414487FF")


int.idents <- names(subsets.all)
i <- 27
png(paste("Figures/Hypo_integrated_subset_", int.idents[[i]], ".png", sep = ""), height = 2.5, width = 7.5, res = 250, units = "in")
plot_grid(DimPlot(object = subsets.all[[i]], reduction = "tsne", group.by = "species", pt.size = .1, cols = cols) + NoAxes() + NoLegend(), DimPlot(object = subsets.all[[i]], reduction = "tsne", group.by = "species.2", pt.size = .1, cols = cols2) + NoAxes() + NoLegend(), DimPlot(subsets.all[[i]], reduction = "tsne", group.by = "integrated_snn_res.0.2", pt.size = .1, label = TRUE) + NoAxes() + NoLegend(), nrow = 1)
dev.off()



# Save Species, species.2, subtype and subcluster plots
png("Figures/Hypo_integrated_tsne_batch_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "basetsne", group.by = "species", pt.size = .1, label = F, cols = cols) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_species.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "species", pt.size = .1, label = F, cols = cols) + NoAxes()+ theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_species2.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "species.2", pt.size = .1, label = F, cols = cols2) + NoAxes() + theme(legend.position = c(0.75,0.95))
dev.off()

png("Figures/Hypo_integrated_tsne_subtype.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subtype", pt.size = .1, label = TRUE) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_subtype_nolabel.png", height = 7.5, width = 7.5, res = 1000, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subtype", pt.size = .1, label = F) + NoAxes() + NoLegend()
dev.off()

png("Figures/Hypo_integrated_tsne_subclustertype_nolabel.png", height = 7.5, width = 7.5, res = 500, units = "in")
DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_SubclusterType", pt.size = .1, label = F) + NoAxes() + NoLegend()
dev.off()

# Whatever you want to split by, should be a factor (not a character vector)
png("Figures/Hypo_integrated_dotplot_subtypes.png", height = 10, width = 9, res = 250, units = "in")
DotPlot(hypo.integrated, features = c("hbaa1", "hbaa2", "hopx", "epd", "gng3", "gad1b", "slc17a6a", "cd74a", "mrc1a", "dusp2", "apoeb", "traf4a", "mpz", "prdx1", "her6"), group.by = "integrated_Subtype", split.by = "species.2", cols = c("blue", "red")) + RotatedAxis()
dev.off()


# Plot groups of genes based on pattern
DotPlot(hypo.integrated, features = rownames(hypo.integrated@assays$RNA[grep("sox", rownames(hypo.integrated@assays$RNA))]), group.by = "integrated_Subtype", split.by = "species.2", cols = c("blue", "red")) + RotatedAxis()

DotPlot(subset(hypo.integrated, idents = "Progenitors"), features = rownames(hypo.integrated@assays$RNA[grep("sox", rownames(hypo.integrated@assays$RNA))]), group.by = "integrated_SubclusterType", split.by = "species.2", cols = c("blue", "red")) + RotatedAxis()

Idents(hypo.integrated) <- "int.Subtype"
DotPlot(subset(hypo.integrated, idents = "Neuronal"), features = c("gng3", "slc17a6a", "gad1b", "rtn4rl2a", "cd9b", "meis2a", "cxcl14", "jdp2b", "penkb", "ngb", "prelid3b", "hspb1", "hspa4a"), group.by = "integrated_Subtype") + RotatedAxis() + coord_flip()


p <- FeaturePlot(object = hypo.integrated, features = c("gng3", "slc17a6b", "gad2", "her6", "prdx1", "otpa", "pfn1", "mpz", "mrc1a", "epd", "hopx", "hbaa1"), cols = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, order = T, combine = F)

for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

png("Figures/hypo_integrated_featureplots_Subtype_markers.png", height = 4.5, width = 12, units = "in", res = 250) 
cowplot::plot_grid(plotlist = p, ncol = 6)
dev.off()


hypo <- hypo.integrated

#### Make tables of cell type proportions

dendrograms <- readRDS(file = "Zeb_Ast_dendrograms_integrated.rds")
orders <- unlist(lapply(dendrograms[[1]], function(x) labels(x)))


prop.table <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species.2)

prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))

prop.table$cell_type <- row.names(prop.table)
prop.table$Subtype <- unique(hypo.integrated@meta.data[, c("integrated_SubclusterType", "integrated_Subtype")])[1:193,2]
prop.table$res.0.2 <- unique(hypo.integrated@meta.data[, c("integrated_SubclusterType", "integrated_SubclusterType_number")])[1:193,2]
# prop.table <- prop.table[,c(5,1,2,3,4)]
prop.table <- melt(prop.table)
prop.table <- prop.table[order(match(prop.table$Subtype, orders)),]

prop.table$cell_type <- factor(prop.table$cell_type, levels = rev(unique(prop.table$cell_type)))


# Plot the dendrogram

offsets <- rev(c(1, 9, 7, 7, 6, 3, 7, 5, 3, 8, 3, 5, 8, 4, 8, 6, 7, 4, 2, 3, 2, 3, 13, 5, 4, 5, 7, 6, 3, 1, 5, 1, 6, 1, 2, 3)*2)
offsets <- offsets + 10

dendro.plot <- ggtree(as.phylo(dendrograms[[1]]), branch.length = "none", layout = "circular") + theme_tree(bgcolor = NA)
dendro.plot <- rotate_tree(dendro.plot, 80) + theme(plot.background = element_blank())
# dendro.plot <- dendro.plot + geom_tiplab2(aes(angle=angle),offset = offsets)


## Do it myself, without the ggtree wrapper - MUCH BETTER

df <- dendro.plot$data
width = 3

prop.table.filt <- prop.table[prop.table$variable == "zebrafish",]
prop.table.filt$y <- df$y[match(prop.table.filt$Subtype, df$label)]
prop.table.filt$col <- ifelse(prop.table.filt$value < 0.1 | prop.table.filt$value > 0.9, 1, 0)
prop.table.filt$size <- ifelse(prop.table.filt$value < 0.1 | prop.table.filt$value > 0.9, 0, 1)
prop.table.filt$x <- max(df$x) + 3 + as.numeric(prop.table.filt$res.0.2) * width

p <- dendro.plot + geom_tile(data = prop.table.filt, aes(x, y, fill = value), width = width, inherit.aes = FALSE) + scale_fill_gradient(low = "#FDE725FF", high = "#414487FF", na.value = NA) + geom_tile(data = prop.table.filt, aes(x, y,color = col, size = size), fill = NA, width = width, inherit.aes = FALSE) + scale_color_gradient(low = NA, high = "red", na.value = "white")

## Add text for cell type labels and for cluster numbers

mapping <- data.frame(label = prop.table.filt$res.0.2, x = prop.table.filt$x)
mapping$y <- prop.table.filt$y

mapping2 <- unique(data.frame(label = prop.table.filt$Subtype, y = prop.table.filt$y))
mapping2$x <- apply(mapping2, 1, function(x) max(prop.table.filt$x[prop.table.filt$Subtype == x[1]]) + 5)
mapping2$angle <- apply(mapping2, 1, function(x) max(dendro.plot$data$angle[dendro.plot$data$label == x[1]], na.rm = T))
mapping2$angle[1:13] <- mapping2$angle[1:13] + 180
mapping2$hjust <- c(rep(1, 13), rep(0, 12))

new.row <- data.frame(label = "Subclusters", y = 27, x = mean(as.numeric(mapping2$x)), angle = 437.5, hjust = .75)

p2 <- p + theme(legend.position = "right", legend.title = element_blank()) + guides(size = FALSE, color = FALSE) + geom_text(data = mapping2, aes(x = x, y = y, label = label), angle = mapping2$angle, color = "black", size = 4, inherit.aes = FALSE, hjust = mapping2$hjust) + geom_text(data = mapping, aes(x = x, y = y, label = label), color = "white", size = 3, inherit.aes = FALSE, hjust = 0.5) 
+ geom_text(data = new.row, aes(x = x, y = y, label = label), angle = new.row$angle, color = "black", size = 6, inherit.aes = FALSE, hjust = new.row$hjust, vjust = .8)



### Make a mosaic plot

ggplot(data = props) + geom_mosaic(aes(x = product(Var1), fill = Var2))


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

pdf("/Volumes/BZ/Home/gizevo30/Manuscripts/Cavefish_single_cell/Figure_4/barplots.pdf", width = 10, height = 3)
ggplot(props, aes(x = Var1, y = value, width = (type_total), fill = Var2)) + geom_bar(stat = "identity", position = "fill", colour = "black", size = 0) + facet_grid(~Var1, scales = "free_x", space = "free_x") + scale_fill_manual(values = cols) + theme(panel.spacing.x = unit(.5, "mm"), strip.text.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 5), axis.title.y = element_text(size = 8), legend.title = element_blank()) + ylab("Prop. of total species cells") + xlab(element_blank()) + geom_hline(yintercept = c(0.3333, 0.6666), linetype = "dashed", color = "black")
dev.off()


hypo.integrated <- ScaleData(object = hypo.integrated, verbose = FALSE)
hypo.integrated <- RunPCA(object = hypo.integrated, npcs = 100, verbose = FALSE)
# hypo.integrated <- RunUMAP(object = hypo.integrated, reduction = "pca", dims = 1:100)
hypo.integrated <- RunTSNE(object = hypo.integrated, reduction = "pca", dims = 1:100, check_duplicates = F)


save(hypo.integr
