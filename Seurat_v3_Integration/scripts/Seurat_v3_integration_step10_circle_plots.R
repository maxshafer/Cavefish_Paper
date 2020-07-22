library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggdendro)
library(phylogram)
library(ggtree)
library(grid)
library(ggmap)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")
hypo <- hypo.integrated
cols0 <- c("lightgoldenrod1", "springgreen4","skyblue2")

# Plot the dendrogram

dendrograms <- readRDS("Zeb_Ast_dendrograms_integrated.rds")
dendrograms <- lapply(dendrograms, function(x) as.dendrogram(x))
 
dendro.plot <- ggtree(as.phylo(dendrograms[[2]]), branch.length = "none", layout = "circular") + theme_tree(bgcolor = NA) # + geom_tiplab(aes(angle=angle))
dendro.plot <- rotate_tree(dendro.plot, 89)
# dendro.plot <- ggplot(as.cladogram(dendrograms[[2]]), labels = FALSE) + coord_polar() + scale_y_reverse(expand = c(.01, 0))

orders <- unlist(lapply(dendrograms[[2]], function(x) labels(x)))


# Plot the # DE genes

species.subclustertype.markers <- readRDS("Hypo_integrated_markers.species.SubclusterType.list.rds")

de.num <- lapply(species.subclustertype.markers, function(x) nrow(x))
de.num[grep("NULL", de.num)] <- 0
de.num <- unlist(de.num)
de.num <- as.data.frame(de.num)

de.num$cell_type <- rownames(de.num)
de.num$z <- 1
de.num$ssc <- "no"
de.num$ssc[grep("species specific cell type", species.subclustertype.markers)] <- "yes"

de.num$cell_type <- factor(de.num$cell_type, levels = rev(orders))


de.circle.plot <- ggplot(de.num, aes(x = cell_type, y = z, fill = de.num, color = ssc)) + geom_tile(aes(fill = de.num)) + scale_fill_gradient2(low = "#4682B4", mid = "#FFFFFF", high = "#FF0000", midpoint = 1, space = "Lab", na.value = "grey50", guide = "colourbar") + scale_color_manual(values = c("transparent", "black")) + theme(text = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + ylim(-8,2) + coord_polar() + guides(fill = FALSE, color = FALSE)


# Proportion figures for sex and orig.ident

ids <- c("integrated_Subtype", "integrated_SubclusterType")

my_colour_palette <- list()
for (i in 1:length(ids)) {
	Idents(hypo.integrated) <- ids[[i]]
	colours <- hue_pal()(length(levels(Idents(hypo))))
	names(colours) <- levels(Idents(hypo))
	my_colour_palette[[i]] <- colours
}
my_colour_palette[[3]] <- my_colour_palette[[1]]
my_colour_palette[[4]] <- my_colour_palette[[2]]

# Make tables of cell type proportions

prop.table <- table(hypo@meta.data$integrated_SubclusterType, hypo@meta.data$species)


prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))

prop.table$Morph_specific <- ifelse(prop.table$zebrafish > .95 | prop.table$zebrafish < 0.05, "yes", "no")

prop.table$cell_type <- row.names(prop.table)
# prop.table <- prop.table[,c(5,1,2,3,4)]
prop.table <- melt(prop.table)
prop.table$cell_type <- as.factor(prop.table$cell_type)


prop.table$cell_type <- factor(prop.table$cell_type, levels = rev(orders))

## Plot Circular barplots!

label_data <- prop.table %>% group_by(cell_type) %>% summarize(tot = sum(value))

label_data$x <- as.numeric(row.names(label_data))
angle = 90-360 * (label_data$x-0.5) / nrow(label_data)
label_data$hjust <- ifelse(angle < -90, 1, 0)
label_data$angle <- ifelse(angle < -90, angle+180, angle)

# Plot the circle dendrogram, with geom_text to add labels
# For SubclusterTypes

prop.circle.plot <- ggplot(prop.table, aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity", aes(color = Morph_specific), size = 1) + theme(text = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + ylim(-3,1) + coord_polar() + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = cols0) + ylab("Species morph Subcluster frequency") + guides(fill = FALSE, color = FALSE) + geom_text(data = label_data, aes(x = as.factor(cell_type), y = tot, label = cell_type, hjust = hjust), colour = "black", size = 4, angle = label_data$angle, inherit.aes = FALSE) + scale_colour_manual(values = c("transparent", "black"))

p <- plot_grid(prop.circle.plot)

p <- p + ggmap::inset(ggplotGrob(de.circle.plot), xmin = 0.105, xmax = .895, ymin = 0.105, ymax = .895)

# dendro.plot + ggmap::inset(ggplotGrob(prop.circle.plot), xmin =0, xmax = 100, ymin = 0, ymax = 100)

pdf("Figures/Hypo_integrated_circle_test.pdf", height = 15, width = 15)
p + ggmap::inset(ggplotGrob(dendro.plot), xmin = 0, xmax = 1, ymin = 0.17, ymax = .83)
dev.off()
