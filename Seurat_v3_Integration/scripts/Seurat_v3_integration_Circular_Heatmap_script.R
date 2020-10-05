library(Seurat)
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
library(dendextend)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

# load objects

hypo <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v3.rds")

cols3 <- c("#FDE725FF", "#22A884FF", "#414487FF")

# Plot the dendrogram

dendrograms <- readRDS("Zeb_Ast_dendrograms_integrated.rds")
#dendrograms <- lapply(dendrograms, function(x) as.dendrogram(x))

# offsets <- rev(c(1, 9, 7, 7, 6, 3, 7, 5, 3, 8, 3, 5, 8, 4, 8, 6, 7, 4, 2, 3, 2, 3, 13, 5, 4, 5, 7, 6, 3, 1, 5, 1, 6, 1, 2, 3)*2)
# offsets <- offsets + 10

dendro.plot <- ggtree(as.phylo(as.dendrogram(dendrograms[[1]])), branch.length = "none", layout = "circular") + theme_tree(bgcolor = NA)
# dendro.plot <- rotate_tree(dendro.plot, 80) + theme(plot.background = element_blank())
# dendro.plot <- dendro.plot + geom_tiplab2(aes(angle=angle),offset = offsets)

orders <- unlist(lapply(dendrograms[[1]], function(x) labels(x)))

# Proportion figures for sex and orig.ident
# Make tables of cell type proportions
# Need to make matrix, for Subtype x Subcluster # (and fill in blanks, in reverse)

# 3 colours
prop.table <- table(hypo@meta.data$integrated_SubclusterType, hypo@meta.data$species.2)
prop.table <- prop.table/colSums(prop.table)

prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))

prop.table$cell_type <- row.names(prop.table)
index <- unique(hypo@meta.data[, c("integrated_SubclusterType", "integrated_Subtype")])
prop.table$Subtype <- index[match(prop.table$cell_type, index[,1]), 2]
prop.table <- prop.table[!(is.na(prop.table$Subtype)),]

index <- unique(hypo@meta.data[, c("integrated_SubclusterType", "integrated_SubclusterType_number")])
prop.table$res.0.4 <- index[match(prop.table$cell_type, index[,1]), 2]

# prop.table$colour <- rgb(1-prop.table$astyanax_cave, 1-prop.table$astyanax_surface, 0, 1)

# prop.table <- prop.table[,c(5,1,2,3,4)]
prop.table <- melt(prop.table)
prop.table <- prop.table[order(match(prop.table$Subtype, orders)),]
prop.table$cell_type <- factor(prop.table$cell_type, levels = rev(unique(prop.table$cell_type)))

## Do it myself, without the ggtree wrapper - MUCH BETTER

df <- dendro.plot$data
width = 3

prop.table.filt <- prop.table[prop.table$variable == "zebrafish",]
prop.table.filt$y <- df$y[match(prop.table.filt$Subtype, df$label)]
prop.table.filt$col <- ifelse(prop.table.filt$value < 0.1 | prop.table.filt$value > 0.9, 1, 0)
prop.table.filt$size <- ifelse(prop.table.filt$value < 0.1 | prop.table.filt$value > 0.9, 0, 1)
prop.table.filt$x <- max(df$x) + 3 + as.numeric(prop.table.filt$res.0.4) * width

# If you want to do a 2 colour graph
p.2 <- dendro.plot + geom_tile(data = prop.table.filt, aes(x, y, fill = value), width = width, inherit.aes = FALSE)  + scale_fill_viridis(option = "D", direction = -1) + geom_tile(data = prop.table.filt, aes(x, y,color = col, size = size), fill = NA, width = width, inherit.aes = FALSE) + scale_color_gradient(low = NA, high = "red", na.value = "white")


## Add text for cell type labels and for cluster numbers

mapping <- data.frame(label = prop.table.filt$res.0.4, x = prop.table.filt$x)
mapping$y <- prop.table.filt$y
mapping2 <- unique(data.frame(label = prop.table.filt$Subtype, y = prop.table.filt$y))
mapping2$x <- apply(mapping2, 1, function(y) max(prop.table.filt$x[prop.table.filt$Subtype == y[1]]) + 5)
mapping2$hjust <- 0
mapping2$hjust[mapping2$y %in% c(1:13)] <- 1
mapping2$angle <- apply(mapping2, 1, function(x) max(dendro.plot$data$angle[dendro.plot$data$label == x[1]], na.rm = T))
mapping3 <- mapping2[order(mapping2$angle),]
mapping3$angle <- sapply(seq_along(c(1:nrow(mapping3))), function(x) min(mapping3$angle) + (x-1)*(max(mapping3$angle) - min(mapping3$angle))/nrow(mapping3))
mapping2 <- mapping3[row.names(mapping2),]

mapping2$angle <- mapping2$angle + 80 # so I can rotate the tree later

mapping2$angle[mapping2$y %in% c(1:13)] <- mapping2$angle[mapping2$y %in% c(1:13)] + 180

new.row <- data.frame(label = "Subclusters", y = 26, x = mean(as.numeric(mapping2$x)), angle = 80, hjust = .75)

# Plot
p2 <- p.2 + theme(legend.position = "right", legend.title = element_blank()) + guides(size = FALSE, color = FALSE) + geom_text(data = mapping2, aes(x = x, y = y, label = label), angle = mapping2$angle, color = "black", size = 6/2.856, inherit.aes = FALSE, hjust = mapping2$hjust) + geom_text(data = mapping, aes(x = x, y = y, label = label), color = "white", size = 6/2.856, inherit.aes = FALSE, hjust = 0.5) + geom_text(data = new.row, aes(x = x, y = y, label = label), angle = new.row$angle, color = "black", size = 10/2.856, inherit.aes = FALSE, hjust = new.row$hjust, vjust = .8)
p2 <- rotate_tree(p2, 80) + theme(plot.background = element_blank())

p2 <- p2 + plot_layout(width = unit(c(135), c("mm")), height = unit(c(135), c("mm")))

dev.new()
p2

