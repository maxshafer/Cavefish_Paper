library(Matrix)
library(dplyr)
library(phylogram)
library(ggtree)
library(ape)
library(cowplot)
library(ggplot2)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k.rds")


#### DF for # of DE genes
#### I don't plot this, and also, these are the species-morph specific marker genes (not DE genes)

gene.lists.pos <- readRDS("marker_gene_lists_pos.rds")

gene.lists.pos.cave <- lapply(seq_along(gene.lists.pos[[1]]), function(x) setdiff(row.names(gene.lists.pos[[3]][[x]]), row.names(gene.lists.pos[[1]][[x]])))
gene.lists.pos.surface <- lapply(seq_along(gene.lists.pos[[1]]), function(x) setdiff(row.names(gene.lists.pos[[2]][[x]]), row.names(gene.lists.pos[[1]][[x]])))

gene.lists.pos.cave.sub <- lapply(seq_along(gene.lists.pos[[4]]), function(x) setdiff(row.names(gene.lists.pos[[6]][[x]]), row.names(gene.lists.pos[[4]][[x]])))
gene.lists.pos.surface.sub <- lapply(seq_along(gene.lists.pos[[4]]), function(x) setdiff(row.names(gene.lists.pos[[5]][[x]]), row.names(gene.lists.pos[[4]][[x]])))


#### Proportion tables for morph and cave

prop.table.morph <- table(hypo@meta.data$SubclusterType, hypo@meta.data$species)
prop.table.morph <- as.data.frame(t(apply(prop.table.morph, 1, function(y) {y/sum(y)})))
prop.table.morph$cell_type <- row.names(prop.table.morph)
prop.table.morph <- reshape2::melt(prop.table.morph)
prop.table.morph$cell_type <- as.factor(prop.table.morph$cell_type)
prop.table.morph$cell_type <- factor(prop.table.morph$cell_type, levels = levels(hypo@meta.data$SubclusterType))

prop.table.cave <- table(hypo@meta.data$SubclusterType, hypo@meta.data$morph)
prop.table.cave <- prop.table.cave[,2:4]
prop.table.cave <- as.data.frame(t(apply(prop.table.cave, 1, function(y) {y/sum(y)})))
prop.table.cave$cell_type <- row.names(prop.table.cave)
prop.table.cave <- reshape2::melt(prop.table.cave)
prop.table.cave$cell_type <- as.factor(prop.table.cave$cell_type)
prop.table.cave$cell_type <- factor(prop.table.cave$cell_type, levels = levels(hypo@meta.data$SubclusterType))


#### DF for the Dissimilarity from the dendrogram - dendrogram distance

dendrograms <- readRDS("Ast_dendrograms.rds")

# phylo <- cophenetic.phylo(as.phylo(dendrograms[["species_SubclusterType"]]))

phylo <- as.data.frame(cophenetic.phylo(as.phylo(dendrograms[[4]]))) # this is missing species-specific subclusters
phylo <- as.matrix(phylo)

# find overlap between cave and surface entries in phylo - this is what I want to extract

test1 <- row.names(phylo)[grep("surface", row.names(phylo))]
test2 <- colnames(phylo)[grep("cave", colnames(phylo))]
test1.sub <- as.character(hypo@meta.data$SubclusterType[match(test1, hypo@meta.data$species_SubclusterType)])
test2.sub <- as.character(hypo@meta.data$SubclusterType[match(test2, hypo@meta.data$species_SubclusterType)])
test3 <- intersect(test1.sub, test2.sub)
test1.index <- match(test3, test1.sub)
test2.index <- match(test3, test2.sub)
test1 <- test1[test1.index]
test2 <- test2[test2.index]

phylo <- phylo[test1, test2]
phylo <- as.matrix(phylo)
phylo.dist <- data.frame(cell_type = test1.sub[test1.index], phylo.dist = diag(phylo))
phylo.dist <- rbind(phylo.dist, data.frame(cell_type = levels(hypo@meta.data$SubclusterType)[!(levels(hypo@meta.data$SubclusterType) %in% phylo.dist$cell_type)], phylo.dist = 0))
phylo.dist$cell_type <- factor(phylo.dist$cell_type, levels = levels(hypo@meta.data$SubclusterType))
phylo.dist$group2 <- "all"

#### Drift Index

DI.sub <- readRDS(file = "Ast_DI.sub.rds")

DI.sub$group <- "all"
df <- data.frame(SubclusterType = levels(hypo@meta.data$SubclusterType)[!(levels(hypo@meta.data$SubclusterType) %in% DI.sub$SubclusterType)], values = NA, group = "all")
df$Subtype <- hypo@meta.data$Subtype[match(df$SubclusterType, hypo@meta.data$SubclusterType)]
DI.sub <- rbind(DI.sub, df)
colnames(DI.sub) <- c("cell_type", "DI", "group", "group2")

## Need to change the x factor levels for combined.2 and label_data
DI.sub$cell_type <- factor(DI.sub$cell_type, levels = levels(hypo@meta.data$SubclusterType))









#### Weir genes plot

weir.genes <- readRDS(file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/astyanax_variants/weir.genes.rds")
snp.genes <- Reduce(union, lapply(weir.genes, function(x) x$gene_name))

combined.markers <- lapply(names(gene.lists.pos$conserved.markers.sub), function(x) unlist(Reduce(union, list(row.names(gene.lists.pos$conserved.markers.sub[[x]]), row.names(gene.lists.pos$cave.markers.sub[[x]]), row.names(gene.lists.pos$surface.markers.sub[[x]])))))
snp.markers <- lapply(seq_along(gene.lists.pos$conserved.markers.sub), function(x) snp.genes[snp.genes %in% combined.markers[[x]]])
per.snp.markers <- unlist(lapply(snp.markers, function(x) length(x)))/unlist(lapply(combined.markers, function(x) length(x)))
per.snp.markers <- data.frame(cell_type = names(gene.lists.pos$conserved.markers.sub), per.snp = per.snp.markers*100, group = hypo@meta.data$Subtype[match(names(gene.lists.pos$conserved.markers.sub), hypo@meta.data$SubclusterType)], group2 = "all")

# Change factor levels
per.snp.markers$group <- factor(per.snp.markers$group, levels = levels(hypo@meta.data$Subtype))
per.snp.markers$cell_type <- factor(per.snp.markers$cell_type, levels = levels(hypo@meta.data$SubclusterType))







#### Label circle plot

label_data <- data.frame(cell_type = levels(hypo@meta.data$SubclusterType), tot = 1, x = seq(1:length(levels(hypo@meta.data$SubclusterType))), angle = 90, hjust = 1, group = hypo@meta.data$Subtype[match(levels(hypo@meta.data$SubclusterType), hypo@meta.data$SubclusterType)])
label_data$cell_type <- factor(label_data$cell_type, levels = levels(prop.table.morph$cell_type))

label_data$x2 <- apply(label_data, 1, function(x) mean(label_data$x[label_data$group == x[6]]))
label_data$subcluster.num <- unlist(lapply(levels(label_data$group), function(x) seq(1:length(label_data$group[label_data$group == x]))))-1

# Add rows to make room for legend, and calculate angles
label_data <- rbind(label_data, data.frame(cell_type = paste("Label", c(1:10), sep = "_"), x = c(168:177), x2 = 176, tot = 0, hjust = 1, angle = 90, subcluster.num = "", group = "Label"))
label_data$angle <- 90-360 * (label_data$x-0.5) / nrow(label_data)
label_data$hjust <- ifelse(label_data$angle < -90, 1, 0)
label_data$angle2 <- apply(label_data, 1, function(x) mean(label_data$angle[label_data$group == x[6]]))
label_data$angle2 <- ifelse(label_data$angle2 < -90, label_data$angle2+180, label_data$angle2)

## Combined all data into one df

# DI.sub$group <- hypo@meta.data$Subtype[match(DI.sub$cell_type, hypo@meta.data$SubclusterType)]

combined <- merge(merge(merge(label_data, DI.sub, sort = F), phylo.dist, sort = F), per.snp.markers, sort = F)
#combined$phylo.dist <- combined$phylo.dist
combined$cell_type <- factor(combined$cell_type, levels = levels(label_data$cell_type))
combined$x <- as.numeric(combined$cell_type)

# The issue is that for Subtypes with only 1 subcluster, it cannot make an area under the plot - can we pad them for x +0.5, x -0.5?
# hmmm, dunno if that will work, or look wonky
# could find the first and last entry, then add another entry x + 0.5, and x - 0.5 that is the same

new.rows <- list()
for (i in 1:length(levels(combined$group))) {
	rows <- (grep(levels(combined$group)[i], combined$group))
	new.row.1 <- combined[min(rows),]
	new.row.1$x <- min(rows - 0.5)
	new.row.2 <- combined[max(rows),]
	new.row.2$x <- max(rows + 0.5)
	
	new.rows[[i]] <- rbind(new.row.1, new.row.2)
}

combined <- rbind(combined, Reduce(rbind, new.rows[1:35]))


#### Add rows to make space for legend

prop.table.morph <- rbind(prop.table.morph, data.frame(cell_type = paste("Label", c(1:10), sep = "_"), variable = c(rep("astyanax_cave", 10), rep("astyanax_surface", 10)), value = 0))
prop.table.cave <- rbind(prop.table.cave, data.frame(cell_type = paste("Label", c(1:10), sep = "_"), variable = c(rep("Molino_cave", 10), rep("Tinaja_cave", 10), rep("Pachon_cave", 10)), value = 0))

DI.sub <- rbind(DI.sub, data.frame(cell_type = paste("Label", x = c(1:10), sep = "_"), DI = NA, group = "Label", group2 = "all"))
per.snp.markers <- rbind(per.snp.markers, data.frame(cell_type = paste("Label", x = c(1:10), sep = "_"), per.snp = NA, group = "Label", group2 = "all"))

combined <- rbind(combined, data.frame(cell_type = paste("Label", x = c(1:10), sep = "_"), group2 = "all", group = "Label", tot = 0.1, x = c(168:177), hjust = 1, angle = 90, x2 = 173, angle2 = 90, subcluster.num = "", DI = NA,phylo.dist = NA, per.snp = NA))




# Plot the circle plots
# For SubclusterTypes

# line.plots <- ggplot(combined, aes(x = x, group = group)) + geom_ribbon(aes(ymin = 0, ymax = de.num, fill = group), alpha = 0.3) + geom_line(data = combined.2, aes(x = x, y = de.num, group = group2), colour = "red") + geom_ribbon(aes(ymin = 0, ymax = phylo.dist, fill = group), alpha = 0.3) + geom_line(data = combined.2, aes(x = x, y = phylo.dist, group = group2), colour = "black") + coord_polar() + ylim(-5000,1615) + theme(text = element_blank(), line = element_blank(), panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA)) + guides(fill = F, colour = F) + geom_hline(yintercept = 500, color = "grey75", size = 1, linetype = "dashed")

phylo.plot <- ggplot(combined, aes(x = x, group = group)) + ylim(-9000,2500) + geom_hline(yintercept = c(750,1500), color = "grey75", size = 1, linetype = "dashed") + geom_ribbon(aes(ymin = 0, ymax = phylo.dist, fill = group), alpha = 0.3) + geom_line(data = combined, aes(x = x, y = phylo.dist, group = group2)) + coord_polar() + theme(text = element_blank(), line = element_blank(), panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA)) + guides(fill = F, colour = F)

# de.num.plot <- ggplot(combined, aes(x = x, group = group, fill = group)) + ylim(-5500,2000) + geom_hline(yintercept = c(500, 1000), color = "grey75", size = 1, linetype = "dashed") + geom_ribbon(aes(ymin = 0, ymax = de.num, fill = group), alpha = 0.3) + geom_line(data = combined.2, aes(x = x, y = de.num, group = group2), colour = "red") + coord_polar() + theme(text = element_blank(), line = element_blank(), panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA)) + guides(fill = F, colour = F)

label.plot <- ggplot(label_data, aes(y=tot, x= cell_type, fill = group, colour = group)) + geom_bar(stat = "identity") + ylim(-1.25,4) + geom_text(data = label_data, aes(x = x2, y = 0.5, label = group), colour = "black", size = 2, hjust = 0.5, angle = label_data$angle2, inherit.aes = FALSE) + theme(text = element_blank(),line = element_blank(), panel.background = element_rect(fill = "transparent", color = NA), plot.background = element_rect(fill = "transparent", color = NA)) + coord_polar() + guides(fill = F, colour = F) + geom_text(data = label_data, aes(x = x, y = 1.1, label = subcluster.num), colour = "black", size = 2, hjust = 0.5, angle = label_data$angle2, inherit.aes = FALSE)

DI.plot <- ggplot(DI.sub, aes(y = DI, x = cell_type, group = group)) + geom_point(colour = "blue", size = 2) + geom_line(colour = "blue", size = 1.5) + coord_polar() + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + theme(text = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA)) + ylim(-5,2) + guides(fill = F)

snp.plot <- ggplot(per.snp.markers, aes(y = .25, x = cell_type, fill = per.snp, group = group)) + geom_tile() + coord_polar() + theme(axis.text = element_blank(), axis.title = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA)) + scale_fill_gradient(low = "white", high = "orange") + ylim(-20,6) + guides(fill = F)

prop.circle.plot <- ggplot(prop.table.morph, aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity") + theme(axis.text = element_blank(), axis.title = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA)) + ylim(-3.25,2) + coord_polar() + geom_hline(yintercept = 0.5, color = "white", size = 1, linetype = "dashed") + ylab("Species morph Subcluster frequency") + scale_fill_viridis_d(direction = -1) + guides(fill = F)

prop.circle.plot.caves <- ggplot(prop.table.cave, aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity") + theme(axis.text = element_blank(), axis.title = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA)) + ylim(-5.5,2) + coord_polar() + geom_hline(yintercept = c(0.333, 0.666), color = "black", size = 0.5, linetype = "dashed") + scale_fill_manual(values = c("#66c2a5", "#fc8d62", "#8da0cb")) + guides(fill = F)


# ## Patchwork version
# 
# # t and b, top and bottom bounds of the area in the grid
# # l and r, left and right bounds
# 
# 
# layout <- c(
#   area(t = 30, l = 30, b = 90, r = 90), # label.plot
#   area(t = 32, l = 32, b = 88, r = 88), # prop.circle.plot
#   area(t = 28, l = 28, b = 92, r = 92), # prop.circle.plot.caves
#   area(t = 22, l = 22, b = 98, r = 98), # DI.plot
#   area(t = 20, l = 20, b = 100, r = 100), # snp.plot
#   area(t = 10, l = 10, b = 110, r = 110)  # phylo.plot
# )
# dev.new()
# label.plot + prop.circle.plot + prop.circle.plot.caves + DI.plot + snp.plot + phylo.plot + plot_layout(design = layout, guides = "collect")


## Use annotation_custom to overlay all of the circle plots

# Need to adjust the snp plot so that it is thinner (change y lims to negative values)

pdf("Figures/hypo_circle_plot_big.pdf", height = 8, width = 8)
plot_grid(label.plot) + annotation_custom(ggplotGrob(prop.circle.plot), xmin = 0.125, xmax = 0.875, ymin = 0.125, ymax = 0.875) + annotation_custom(ggplotGrob(prop.circle.plot.caves), xmin = 0.0625, xmax = 0.9375, ymin = 0.0625, ymax = 0.9375) + annotation_custom(ggplotGrob(DI.plot), xmin = -0.015, xmax = 1.015, ymin = -0.015, ymax = 1.015) + annotation_custom(ggplotGrob(snp.plot), xmin = -0.08, xmax = 1.08, ymin = -0.08, ymax = 1.08) + annotation_custom(ggplotGrob(phylo.plot), xmin = -0.115, xmax = 1.115, ymin = -0.115, ymax = 1.115)
dev.off()




### OLD

# Plot the dendrogram

dendrograms <- hypo@cluster.tree
dendrograms <- lapply(dendrograms, function(x) as.dendrogram(x))
dendro.plot <- ggtree(as.phylo(dendrograms[[2]]), branch.length = "none", layout = "circular") + theme_tree(bgcolor = NA) # + geom_tiplab(aes(angle=angle))
dendro.plot <- rotate_tree(dendro.plot, 89)

## other old

phylo.plot <- ggplot(phylo.dist, aes(y = phylo.dist, x = cell_type,)) + geom_bar(stat = "identity", fill = "black") + coord_polar() + theme(text = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent")) + ylim(-1666, 2500) + guides(fill = F)

phylo.plot <- ggplot(phylo.dist, aes(y = phylo.dist, x = cell_type, group = group)) + geom_point() + geom_line() + coord_polar() + theme(text = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent")) + ylim(-1666, 2500) + guides(fill = F)


de.circle.plot <- ggplot(de.num, aes(x = cell_type, y = de.num, fill = de.num)) + geom_bar(stat = "identity") + coord_polar() + theme(text = element_blank(), line = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_rect(fill = "transparent", color = NA)) + ylim(-5000,1615) + guides(fill = F)


plot_grid(label.plot) + annotation_custom(ggplotGrob(prop.circle.plot), xmin = 0, xmax = 1, ymin = 0, ymax = 1) + annotation_custom(ggplotGrob(prop.circle.plot.caves), xmin = -0.045, xmax = 1.045, ymin = -0.045, ymax = 1.045) + annotation_custom(ggplotGrob(DI.plot), xmin = -0.2, xmax = 1.2, ymin = -0.2, ymax = 1.2) + annotation_custom(ggplotGrob(phylo.plot), xmin = -0.2, xmax = 1.2, ymin = -0.2, ymax = 1.2)


plot_grid(line.plots) + annotation_custom(ggplotGrob(label.plot), xmin = 0, xmax = 1, ymin = 0, ymax = 1)
							+ annotation_custom(ggplotGrob(prop.circle.plot.caves), xmin = 0, xmax = 1, ymin = 0, ymax = 1) 
							+ annotation_custom(ggplotGrob(DI.plot), xmin = -0.045, xmax = 1.045, ymin = -0.045, ymax = 1.045) 
							+ annotation_custom(ggplotGrob(phylo.plot), xmin = 0.025, xmax = 0.975, ymin = 0.025, ymax = 0.975) 



plot_grid(prop.circle.plot) + ggmap::inset(ggplotGrob(prop.circle.plot.caves), xmin = 0, xmax = 1, ymin = 0, ymax = 1) + ggmap::inset(ggplotGrob(de.circle.plot), xmin = -0.045, xmax = 1.045, ymin = -0.045, ymax = 1.045) + ggmap::inset(ggplotGrob(DI.plot), xmin = -0.045, xmax = 1.045, ymin = -0.045, ymax = 1.045)



+ geom_text(data = label_data, aes(x = as.factor(cell_type), y = tot + 1.5, label = cell_type, hjust = hjust), colour = "black", size = 4, angle = label_data$angle, inherit.aes = FALSE)
