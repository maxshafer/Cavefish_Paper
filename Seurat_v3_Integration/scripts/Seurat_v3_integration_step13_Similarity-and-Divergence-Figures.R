library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(viridis)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v3.rds")
gene.lists.pos <- readRDS("drift_gene_lists_pos.rds")

SI.list <- readRDS("SI_results_trinarized_markers.rds")
SI <- SI.list[[1]]
SI.sub <- SI.list[[2]]
SI.ast <- SI.list[[3]]
SI.ast.sub <- SI.list[[4]]

## Make Correlation matrix

norm.cluster <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Normed_expression_data.rds")
str(norm.cluster, max.level = 2)

intersect.genes <- intersect(row.names(norm.cluster[[4]][[1]]), row.names(norm.cluster[[4]][[2]]))

zeb.avg <- norm.cluster[[4]][[1]][intersect.genes,]
ast.avg <- norm.cluster[[4]][[2]][intersect.genes,]

zeb.ast <- cor(zeb.avg, ast.avg, method = "pearson") # row x column
# zeb.ast <- zeb.ast/rowMax(zeb.ast) # scale by row, looks very similar
zeb.ast <- reshape2::melt(zeb.ast)

Corr.matrix.scaled <- ggplot(zeb.ast, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 3)) + scale_fill_viridis(direction = 1, limits = c(0,1))
Corr.matrix.scaled <- Corr.matrix.scaled + theme(axis.title = element_blank())

## Make comparison figure SI between zeb-ast and cave-surface

SI.combined <- rbind(SI, SI.ast)
SI.combined$species <- c(rep("Between species", nrow(SI)), rep("Between species-morphs", nrow(SI.ast)))
species.SI <- ggplot(SI.combined, aes(x = species, y = values, colour = species, group = species)) + geom_jitter() + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic() + ylim(c(0,1)) + scale_colour_viridis_d()
species.SI <- species.SI + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")")))

SI.combined.sub <- rbind(SI.sub[,1:3], SI.ast.sub)
SI.combined.sub$species <- c(rep("Between species", nrow(SI.sub)), rep("Between species-morphs", nrow(SI.ast.sub)))
species.SI.sub <- ggplot(SI.combined.sub, aes(x = species, y = values, colour = species, group = species)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic() + ylim(c(0,1)) + scale_colour_viridis_d()
species.SI.sub <- species.SI.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")"))) + NoLegend()

## Plot Figure S4
dev.new()
Corr.matrix.scaled + plot_layout(height = unit(c(90), c("mm")), width = unit(c(90), c("mm")))
dev.new()
species.SI.sub + plot_layout(height = unit(c(30), c("mm")), width = unit(c(15), c("mm")))

## Similarity Index (SI) and % paralog figures

matrix.ast2.zeb.sub <- SI.list[[11]]
names(matrix.ast2.zeb.sub) <- names(gene.lists.pos$conserved.markers.sub)


## Make row-scaled heatmap
matrix.data <- lapply(matrix.ast2.zeb.sub, function(x) unlist(x))
matrix.data <- reshape2::melt(lapply(matrix.data, function(x) (x)/max(x)))
# matrix.data <- reshape2::melt(matrix.data)
matrix.data$L2 <- rep(names(matrix.ast2.zeb.sub[[1]]),151)
matrix.data$L1 <- factor(matrix.data$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
matrix.data$L2 <- factor(matrix.data$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
SI.matrix.scaled <- ggplot(matrix.data, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 3)) + scale_fill_viridis(direction = 1, limits = c(0,1))

SI.matrix.scaled <- SI.matrix.scaled + theme(axis.title = element_blank())


## Plot progenitors versus neurons

SI.sub$Neuronal <- "Neuronal"
SI.sub$Neuronal[grepl("Progenitors", SI.sub$Subcluster)] <- "Progenitors"
SI.sub$Neuronal <- factor(SI.sub$Neuronal, levels = c("Progenitors", "Neuronal"))
neuronal <- ggplot(as.data.frame(SI.sub), aes(x = Neuronal, y = values, group = Neuronal, color = Neuronal)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic()
neuronal <- neuronal + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")")))


SI.compare <- data.frame(cell_type = SI$Cluster, cluster = SI$values, cell_types = aggregate(SI.sub[,2], list(SI.sub$Cluster), mean)[,2])
cluster <- ggplot(reshape::melt(SI.compare), aes(x = variable, y = value, group = cell_type, color = cell_type)) + geom_jitter(size = 0.5) + geom_line(aes(group = cell_type), colour = "grey75") + geom_boxplot(aes(group = variable),outlier.color = NA, fill = "transparent", color = "black") + theme_classic() + theme(legend.position = "none")
cluster <- cluster + scale_x_discrete(labels=c("cluster" = "Cluster", "cell_type" = "Subcluster")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")"))) 

# t.test(SI.sub[SI.sub$Neuronal == "Neuronal", "values"], SI.sub[SI.sub$Neuronal == "Progenitors", "values"])
# t.test(SI.compare$cell_types, SI.compare$cluster, paired = T)

### Paralog plots
## Make SI by Paralog (sum?? or mean??)

para.results <- readRDS("Paralog_results_trinarized_markers.rds")
para.data <- reshape2::melt(para.results[[1]][,13:18])
para.data.sub <- reshape2::melt(para.results[[2]][,13:19])

percent.para <- ggplot(para.data.sub[para.data.sub$variable == c("percent.para.1", "percent.para.2"),], aes(x = variable, y = value, colour = variable, group = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic() + ylim(c(0,18))
percent.para <- percent.para + scale_x_discrete(labels=c("percent.para.1" = "Zebrafish", "percent.para.2" = "Mexican tetra")) + theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) + ylab("% Paralogs")

odds.para <- ggplot(para.data.sub[para.data.sub$variable == c("fisher.odds.1", "fisher.odds.2"),], aes(x = variable, y = value, colour = variable, group = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic()
odds.para <- percent.para + scale_x_discrete(labels=c("fisher.odds.1" = "Zebrafish", "fisher.odds.2" = "Mexican tetra")) + theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) + ylab("Odds Ratio")

## Divergence (dT) plots
## Figure 2d-e, S4c

dT.list <- readRDS(file = "dT_list.rds")

branch.levels <- levels(dT.list[[2]]$branch)
branch.levels <- c(branch.levels, "Danio rerio")

dT.list[[1]]$branch <- factor(dT.list[[1]]$branch, levels = branch.levels)
dT.list[[2]]$branch <- factor(dT.list[[2]]$branch, levels = branch.levels)

plot.zeb <- ggplot(dT.list[[1]], aes(x = dT, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour") + theme_classic()
plot.ast <- ggplot(dT.list[[2]], aes(x = dT, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour") + theme_classic()

plot.zeb <- plot.zeb + guides(color = guide_legend(title = "LCA of paralogs"))
plot.ast <- plot.ast + guides(color = guide_legend(title = "LCA of paralogs"))

## Compare gene pairs which are common to zebrafish and cavefish
## Compare dTs for cavefish and zebrafish, and correlate binarized gene expression patterns

library(plyr)

dTzeb.2 <- plyr::match_df(dT.list[[1]], dT.list[[2]], on = c("gene1", "gene2"))
dTast.2 <- plyr::match_df(dT.list[[2]], dT.list[[1]], on = c("gene1", "gene2"))
dT.combined <- merge(dTast.2, dTzeb.2, by = c("gene1", "gene2"))

# x is ast, y is zeb
all.genes <- ggplot(dT.combined, aes(x = dT.x, y = dT.y, color = branch.x, label = paste(gene1, gene2, sep = "-"))) + geom_point(colour = "black", size = 0.5) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + theme_classic() + ylab("dT zebrafish") + xlab("dT Mexican tetra") + ggtitle("All gene pairs", subtitle = "Pearson correlation = 0.54")
otophysi <- ggplot(dT.combined[dT.combined$branch.x == "Otophysi" | dT.combined$branch.y == "Otophysi",], aes(x = dT.x, y = dT.y)) + geom_point(colour = "#FF61C9", size = 0.5) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + theme_classic() + ylab("dT zebrafish") + xlab("dT Mexican tetra") + ggtitle("LCA Otophysi gene pairs", subtitle = "Pearson correlation = 0.45")


# Put together Figure 2

part1 <- SI.matrix.scaled + (neuronal + cluster + percent.para + odds.para + plot_layout(ncol = 4)) + plot_layout(nrow = 2, heights = c(3,1), guides = "collect")
part1 <- part1 + plot_layout(height = unit(c(90,30), c("mm", "mm")), width = unit(c(90), c("mm")))
dev.new()
part1

part2 <- plot.zeb + plot.ast + all.genes + otophysi + plot_layout(ncol = 1)
part2 <- part2 + plot_layout(height = unit(c(90/4, 90/4, 90/4, 90/4), c("mm", "mm", "mm", "mm")), width = unit(c(90/4), c("mm")))
dev.new()
part2


### Prep figures for Figure S5

## Make orthology metrics figures

orth_conf_metrics <- readRDS(file = "paralogs-orthology-conf-metrics.rds")

confidence <- ggplot(reshape2::melt(orth_conf_metrics[[1]][,2:3]), aes(x = variable, y = value, group = variable, colour = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black")
confidence <- confidence + theme_classic() + ylim(c(0,1))  + ylab("Gene orthology confidence")
confidence <- confidence + theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + NoLegend()

gene.order <- ggplot(reshape2::melt(orth_conf_metrics[[2]][,2:3]), aes(x = variable, y = value, group = variable, colour = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black")
gene.order <- gene.order + theme_classic() + ylim(c(0,100)) + ylab("Gene order score")
gene.order <- gene.order + theme(axis.title.x = element_blank(), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + NoLegend()

# Make synteny % figure

syn.corr <- readRDS(file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/synteny-corrected-percentage.rds")

syn.corr.plot <- ggplot(reshape2::melt(syn.corr), aes(x = variable, y = value*100, colour = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA) + theme_classic()
syn.corr.plot <- syn.corr.plot + theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + ylab("% Synteny-corrected genes") + NoLegend()
syn.corr.plot <- syn.corr.plot + scale_x_discrete(breaks=c("con.markers","zeb.markers", "ast.markers", "zeb.para.markers","ast.para.markers"), labels=c("Conserved \n marker genes", "Zebrafish \n marker genes", "Mexican tetra \n marker genes", "Zebrafish \n paralogs", "Mexican tetra \n paralogs"))


## Make paralog % by drift scatter

para.results[[2]]$Neuronal <- "Non-Neuronal"
para.results[[2]]$Neuronal[grepl("Neuronal", para.results[[2]]$Subcluster)] <- "Neuronal"

r <- cor((para.results[[2]]$percent.para.2 + para.results[[2]]$percent.para.1)/2, para.results[[2]]$driftindex)
para.results.neuronal <- para.results[[2]][para.results[[2]]$Neuronal == "Neuronal",]
r <- cor((para.results.neuronal$percent.para.2 + para.results.neuronal$percent.para.1)/2, para.results.neuronal$driftindex)
para.results.non <- para.results[[2]][para.results[[2]]$Neuronal == "Non-Neuronal",]
r <- cor((para.results.non$percent.para.2 + para.results.non$percent.para.1)/2, para.results.non$driftindex)

## If you group by neuronal/non, there is a clear difference in the relationship
para.drift <- ggplot(para.results[[2]], aes(x = (percent.para.2 + percent.para.1)/2, y = driftindex, color = Cluster, label = Subclusters)) + geom_point(size = 0.5) + guides(color = F) + ylab("Drift Index") + xlab("Mean of % Paralogs per species") # + geom_smooth(method = lm, se = T, color = "black", fill = "grey85")
para.drift <- para.drift + theme_classic() + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + annotate(geom = "text", x = 2.5, y = 0.55, label = paste("Pearson r = ", round(r,digits = 2), sep = ""), color = "black") + NoLegend()


## Make corrected SI heatmap figures

## Row-scaled delta between SI and Corrected-SI
SI.delta <- lapply(seq_along(SI.list[[11]]), function(x) lapply(seq_along(SI.list[[11]][[x]]), function(y) SI.list[[12]][[x]][[y]] - SI.list[[11]][[x]][[y]]))
names(SI.delta) <- names(gene.lists.pos$conserved.markers.sub)

SI.delta <- lapply(SI.delta, function(x) unlist(x))
SI.delta <- reshape2::melt(lapply(SI.delta, function(x) (x)/max(x)))
SI.delta$L2 <- rep(names(SI.list[[12]][[1]]),151)
SI.delta$L1 <- factor(SI.delta$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
SI.delta$L2 <- factor(SI.delta$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))

SI.delta.plot <- ggplot(SI.delta, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 1.5), axis.title = element_blank()) + scale_fill_viridis(direction = 1, option = "C")

## Make jitter plot comparing uncorrected and corrected SI by subcluster

matrix.data <- lapply(SI.list[[11]], function(x) unlist(x))
names(matrix.data) <- names(gene.lists.pos$conserved.markers.sub)
matrix.data <- reshape2::melt(matrix.data)
matrix.data$L2 <- rep(names(matrix.ast2.zeb.sub[[1]]),151)
matrix.data$L1 <- factor(matrix.data$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
matrix.data$L2 <- factor(matrix.data$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))

matrix.data.corr <- lapply(SI.list[[12]], function(x) unlist(x))
names(matrix.data.corr) <- names(gene.lists.pos$conserved.markers.sub)
matrix.data.corr <- reshape2::melt(matrix.data.corr)
matrix.data.corr$L2 <- rep(names(matrix.ast2.zeb.sub[[1]]),151)
matrix.data.corr$L1 <- factor(matrix.data.corr$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
matrix.data.corr$L2 <- factor(matrix.data.corr$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))

SI.combined.corr <- rbind(matrix.data[matrix.data$L1 == matrix.data$L2,], matrix.data.corr[matrix.data.corr$L1 == matrix.data.corr$L2,])
SI.combined.corr$sample <- c(rep("Uncorrected SI", 151), rep("Corrected SI", 151))
SI.combined.corr$sample <- factor(SI.combined.corr$sample, levels = c("Uncorrected SI", "Corrected SI"))

corr.delta.plot <- ggplot(SI.combined.corr, aes(x = sample, y = value, color = sample)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d(option = "C") + scale_color_manual(values = c("#5D01A6FF", "#FDB32FFF")) + theme_classic()
corr.delta.plot <- corr.delta.plot + theme(axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab("Δ SI vs Corrected SI") + NoLegend()

# Put them together into S5

dev.new()
confidence + gene.order + corr.delta.plot + plot_layout(nrow = 1, widths = unit(c(10,10,10), "mm"), height = unit(c(20,20,20), "mm"), guides = "collect")
dev.new()
syn.corr.plot + para.drift + plot_layout(nrow = 1, widths = unit(c(15,25), "mm"), height = unit(c(25,25), "mm"), guides = "collect")
dev.new()
SI.delta.plot + plot_layout(widths = unit(c(90), "mm"), height = unit(c(90), "mm"), guides = "collect")


## Make Paralog divergence supplemental figure S7

##### Correlate binarized cell types between species for each pair

# Create matrices of normalized expression for each gene across integrated_Subclusters

zeb.int <- Reduce(cbind, norm.cluster[[4]][[1]])
colnames(zeb.int) <- names(norm.cluster[[4]][[1]])
row.names(zeb.int) <- row.names(norm.cluster[[4]][[1]][[1]])

ast.int <- Reduce(cbind, norm.cluster[[4]][[2]])
colnames(ast.int) <- names(norm.cluster[[4]][[2]])
row.names(ast.int) <- row.names(norm.cluster[[4]][[2]][[1]])

dTcombined <- lapply(unique(dT.combined$branch.x), function(x) dT.combined[dT.combined$branch.x == x,])
names(dTcombined) <- unique(dT.combined$branch.x)

# Match the gene pairs from the dTcombined lists from each .int matrix and cor() them

correlations <- lapply(dTcombined, function(x) apply(x, 1, function(x) cor(unlist(c(zeb.int[as.character(x[1]),intersect(colnames(zeb.int), colnames(ast.int))], zeb.int[as.character(x[2]),intersect(colnames(zeb.int), colnames(ast.int))])), unlist(c(ast.int[as.character(x[1]),intersect(colnames(zeb.int), colnames(ast.int))], ast.int[as.character(x[2]),intersect(colnames(zeb.int), colnames(ast.int))])))))

lapply(correlations, function(x) mean(x))

correlations2 <- data.frame(branch = c(rep(names(correlations)[[1]], length(correlations[[1]])), rep(names(correlations)[[2]], length(correlations[[2]])), rep(names(correlations)[[3]], length(correlations[[3]])), rep(names(correlations)[[4]], length(correlations[[4]])), rep(names(correlations)[[5]], length(correlations[[5]])), rep(names(correlations)[[6]], length(correlations[[6]])), rep(names(correlations)[[7]], length(correlations[[7]])), rep(names(correlations)[[8]], length(correlations[[8]]))), value = unlist(correlations))

correlations2$branch <- factor(correlations2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))

anno_df <- compare_means(value ~ branch, data = correlations2)

# Figure S7a
corr.plot <- ggplot(correlations2, aes(x = branch, y = value)) + geom_jitter(color = "grey65", shape = 21) + geom_boxplot(fill = "transparent", outlier.color = NA) + theme_classic() + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + stat_compare_means()


## Do it for binarized expression (there or not)

zeb.int.b <- ifelse(zeb.int > 0, 1, 0)
ast.int.b <- ifelse(ast.int > 0, 1, 0)

correlations.b <- lapply(dTcombined, function(x) apply(x, 1, function(x) cor(unlist(c(zeb.int.b[as.character(x[1]),intersect(colnames(zeb.int.b), colnames(ast.int.b))], zeb.int.b[as.character(x[2]),intersect(colnames(zeb.int.b), colnames(ast.int.b))])), unlist(c(ast.int.b[as.character(x[1]),intersect(colnames(zeb.int.b), colnames(ast.int.b))], ast.int.b[as.character(x[2]),intersect(colnames(zeb.int.b), colnames(ast.int.b))])))))

lapply(correlations.b , function(x) mean(x))

correlations.b.2 <- data.frame(branch = c(rep(names(correlations)[[1]], length(correlations.b[[1]])), rep(names(correlations)[[2]], length(correlations.b[[2]])), rep(names(correlations)[[3]], length(correlations.b[[3]])), rep(names(correlations)[[4]], length(correlations.b[[4]])), rep(names(correlations)[[5]], length(correlations.b[[5]])), rep(names(correlations)[[6]], length(correlations.b[[6]])), rep(names(correlations)[[7]], length(correlations.b[[7]])), rep(names(correlations)[[8]], length(correlations.b[[8]]))), value = unlist(correlations.b))

correlations.b.2$branch <- factor(correlations.b.2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))

anno_df <- compare_means(value ~ branch, data = correlations.b.2)

# Figure S7b
s7b <- ggplot(correlations.b.2[correlations.b.2$branch != "Opisthokonta",], aes(x = branch, y = value)) + geom_jitter(color = "grey65", shape = 21) + geom_boxplot(fill = "transparent", outlier.color = NA) + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + stat_compare_means()

s7b <- s7b + theme_classic() + stat_compare_means()

### Make DotPlots for divergent gene pairs

# hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v2.rds")

Idents(hypo.integrated) <- "species.2"

dot.zeb <- DotPlot(hypo.integrated.zeb, features = c("etv5b", "etv5a"), group.by = "integrated_Cluster", scale.max = 60) + RotatedAxis() + coord_flip() + scale_color_viridis() + theme(axis.title = element_blank(), axis.text = element_text(size = 8), axis.text.x = element_blank())
dot.ast <- DotPlot(hypo.integrated.ast, features = c("etv5b", "etv5a"), group.by = "integrated_Cluster", scale.max = 60) + RotatedAxis() + coord_flip() + scale_color_viridis() + theme(axis.title = element_blank(), axis.text = element_text(size = 8))

# Figure S4d
dots <- dot.zeb / dot.ast + plot_layout(nrow = 2, guides = "collect")

# Put them together

dev.new()
s7b + plot_layout(width = unit(75, "mm"), height = unit(30, "mm"))

dev.new()
dots + plot_layout(width = unit(90, "mm"), height = unit(15, "mm"))














# ##### others, not used
# 
# ## Plot the number of conserved and species specific marker genes per cluster
# 
# test <- lapply(seq_along(1:length(gene.lists.pos[[7]])), function(x) {
#   zeb <- length(setdiff(rownames(gene.lists.pos[[8]][[x]]), rownames(gene.lists.pos[[7]][[x]])))
#   ast <- length(setdiff(rownames(gene.lists.pos[[9]][[x]]), rownames(gene.lists.pos[[7]][[x]])))
#   con <- length(rownames(gene.lists.pos[[7]][[x]]))
#   return(data.frame(con, ast, zeb)) } )
# 
# test <- Reduce(rbind, test)
# test$cell_type <- names(gene.lists.pos[[7]])
# 
# ggplot(melt(test), aes(y = variable, x = cell_type,color = value)) + geom_point(size = 5) + scale_color_viridis(option = "B") + coord_flip() + theme(axis.text = element_text(size = 6))
# 
# ## Plot the # of conserved and species specific markers
# 
# marker.numbers <- data.frame(subtype = names(unlist(lapply(gene.lists.pos[[1]], function(x) nrow(x)))), conserved = unlist(lapply(gene.lists.pos[[1]], function(x) nrow(x))), zebrafish = unlist(lapply(gene.lists.pos[[2]], function(x) nrow(x))), cavefish = unlist(lapply(gene.lists.pos[[3]][!(names(gene.lists.pos[[3]]) %in% "Cilliated")], function(x) nrow(x))))
# marker.numbers.sub <- data.frame(subcluster = names(unlist(lapply(gene.lists.pos[[7]], function(x) nrow(x)))), subtype = hypo.integrated@meta.data$integrated_Cluster[match(names(unlist(lapply(gene.lists.pos[[7]], function(x) nrow(x)))), hypo.integrated@meta.data$integrated_Subcluster)], conserved = unlist(lapply(gene.lists.pos[[7]], function(x) nrow(x))), zebrafish = unlist(lapply(gene.lists.pos[[8]], function(x) nrow(x))), cavefish = unlist(lapply(gene.lists.pos[[9]], function(x) nrow(x))))
# marker.numbers$subtype <- factor(marker.numbers$subtype, levels = levels(hypo.integrated@meta.data$integrated_Cluster))
# marker.numbers.sub$subtype <- factor(marker.numbers.sub$subtype, levels = levels(hypo.integrated@meta.data$integrated_Cluster))
# 
# ggplot(marker.numbers.sub, aes(x = subtype, y = conserved, fill = subtype)) + geom_boxplot(outlier.color = "transparent") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + guides(fill = F) + coord_flip() + geom_jitter(data = marker.numbers.sub, aes(x = subtype, y = conserved), color = "black") + geom_point(data = marker.numbers, aes(x = subtype, y = conserved), color = "red")
# 
# 
# 
# ### Extra plots not used
# # subtypes <- melt(subtypes[,9:12])
# 
# ggplot(subtypes, aes(x = percent.para.2 + percent.para.1, y = driftindex, color = Clusters, label = Clusters)) + geom_point() + geom_text() + guides(color = FALSE) + geom_smooth(method = lm, se = T, color = "black", fill = "grey85")
# ggplot(subtypes, aes(x = percent.para.1, y = percent.para.2, color = Clusters, label = Clusters)) + geom_point() + geom_text() + guides(color = FALSE)
# ggplot(subtypes.2[subtypes.2$variable == "percent.para.1" | subtypes.2$variable == "percent.para.2",], aes(Clusters,x = variable, y = value, color = Clusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F) + ylab("Percentage Paralog") + xlab("Species") + scale_x_discrete(labels = c("D. rerio", "A. mexicanus")) + theme(axis.text.x = element_text(face = "italic"))
# 
# ggplot(subtypes.2[subtypes.2$variable == "fisher.odds.1" | subtypes.2$variable == "fisher.odds.2",], aes(Clusters,x = variable, y = value, color = Clusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F)
# ggplot(subtypes.2[subtypes.2$variable == "driftindex",], aes(Clusters,x = variable, y = value, color = Clusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 3, position = position_jitter(0.1)) + guides(color = F)
# 
# 
# fisher.test(matrix(data = apply(do.call(rbind, paralog.numbers), 2, function(x) sum(x))[5:8], ncol = 2, nrow = 2))
# 
# 
# ## FOR subclusters
# 
# ggplot(subclusters, aes(x = percent.para.2 + percent.para.1, y = driftindex, color = Subclusters, label = Subclusters)) + geom_point() + guides(color = FALSE) + geom_smooth(method = lm, se = T, color = "black", fill = "grey85")
# ggplot(subclusters, aes(x = percent.para.1, y = driftindex, color = Subclusters, label = Subclusters)) + geom_point() + guides(color = FALSE)
# 
# 
# ggplot(subclusters.2[subclusters.2$variable == "percent.para.1" | subclusters.2$variable == "percent.para.2",], aes(Subclusters,x = variable, y = value, color = Subclusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 1, position = position_jitter(0.1)) + guides(color = F)
# ggplot(subclusters.2[subclusters.2$variable == "fisher.odds.1" | subclusters.2$variable == "fisher.odds.2",], aes(Subclusters,x = variable, y = value, color = Subclusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 1, position = position_jitter(0.1)) + guides(color = F)
# ggplot(subclusters.2[subclusters.2$variable == "driftindex",], aes(Subclusters,x = variable, y = value, color = Subclusters, group = variable)) + geom_boxplot(width = 0.5, outlier.shape = NA) + geom_jitter(size = 1, position = position_jitter(0.1)) + guides(color = F)
# 
# 
# fisher.test(matrix(data = apply(do.call(rbind, paralog.numbers.sub), 2, function(x) sum(x))[5:8], ncol = 2, nrow = 2))
# 
# 
# ## Make Patchwork plot for Figure 3
# library(patchwork)
# 
# # Takes plots from the PercentParalogs.R, Drift_Index.R, and expression_divergence_v2.R scripts

( (SI.plot | para.drift) + plot_layout(ncol = 2, widths = c(3,1)) ) / ( (percent.para | plot.id) + plot_layout(ncol = 2, widths = c(3,1)) )
