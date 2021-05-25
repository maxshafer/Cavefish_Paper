library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(viridis)
library(ggpubr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

a = 1.5
b = 2
f = 0.1

## Divergence (dT) plots
## Figure 2d-e, S4c

dT.list <- readRDS(file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

dT.NF <- dT.list[[4]]

# Load 2R and 3R lists and sort pairs for matching

twor <- read.csv("drerio.Pairs.Strict.2R.txt", sep = "\t")
threer <- read.csv("drerio.Pairs.Strict.3R.txt", sep = "\t")

twor <- as.data.frame(t(apply(twor[, c(3,4)], 1, sort)))
threer <- as.data.frame(t(apply(threer[, c(3,4)], 1, sort)))

colnames(twor) <- c("gene1", "gene2")
colnames(threer) <- c("gene1", "gene2")

twor$genes <- paste(twor$gene1, twor$gene2, sep = "_")
threer$genes <- paste(threer$gene1, threer$gene2, sep = "_")

# Add origin as column to dT.NF
dT.NF$genes <- paste(dT.NF$gene1, dT.NF$gene2, sep = "_")

dT.NF$origin <- ifelse(dT.NF$genes %in% twor$genes, "2R", ifelse(dT.NF$genes %in% threer$genes, "3R", "Other"))

# set factor levels for branch and origin, and melt data.frame

branch.levels <- c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Clupeocephala", "Otophysi", "Danio rerio", "Astyanax mexicanus")

dT.NF$branch <- factor(dT.NF$branch, levels = branch.levels)
dT.NF$origin <- factor(dT.NF$origin, levels = c("Other", "3R", "2R"))

# Remove the divide by # of cell types, otherwise the results are too weird, 161 in zeb, 184 in ast
dT.NF$NF.gene1.ast <- dT.NF$NF.gene1.ast*184
dT.NF$NF.gene2.ast <- dT.NF$NF.gene2.ast*184
dT.NF$NF.gene1.zeb <- dT.NF$NF.gene1.zeb*161
dT.NF$NF.gene2.zeb <- dT.NF$NF.gene2.zeb*161

dT.NF.melt <- reshape2::melt(dT.NF)

########################################################################################
########################################################################################
########################################################################################

## Plot ecdf plots for dT, redudancy, and specialiation scores grouped by age (branch)
dt.ecdf <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("dT.zeb", "dT.ast"),], aes(x = value, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour")
dt.ecdf <- dt.ecdf + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + theme_classic() + guides(color = guide_legend(title = "LCA of paralogs")) + xlab("dT")

redun.ecdf <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("dT.zeb.ancestral", "dT.ast.ancestral"),], aes(x = value, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour")
redun.ecdf <- redun.ecdf + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + theme_classic() + guides(color = guide_legend(title = "LCA of paralogs")) + xlab("Redundancy score") + ylab("")

sc.ecdf <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("NF.gene1.zeb", "NF.gene2.zeb", "NF.gene1.ast", "NF.gene2.ast"),], aes(x = value, color = branch)) + stat_ecdf(size = 1) + ylab("% of paralog pairs") + scale_color_brewer(type = "seq", palette = "YlOrRd", direction = 1, aesthetics = "colour")
sc.ecdf <- sc.ecdf + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), legend.text = element_text(size = 8), legend.title = element_text(size = 10)) + theme_classic() + guides(color = guide_legend(title = "LCA of paralogs")) + xlab("Specialization score") + ylab("")

## Plot dT, redundancy, and specialiation scores grouped by origin (2R vs 3R)

# Density plot of dT values, grouped by strict 2R or 3R
# Jitter plot of dT values, grouped by strict 2R or 3R
# 2R vs 3R is significantly different, Wilcoxon p-value = 0.0076
dT.origin.density <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("dT.zeb", "dT.ast") & dT.NF.melt$origin != "Other",], aes(x = value, group = origin, fill = origin, colour = origin)) + geom_density(alpha = 0.1, size = 1) + theme_classic() + scale_colour_viridis_d() + scale_fill_viridis_d() + theme(legend.position = c(0.65,0.75), legend.title = element_blank())
dT.origin.density <- dT.origin.density + ylab("Density") + xlab("log(Ancestral dT)") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), line = element_line(size = 1))

dT.origin.jitter <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("dT.zeb", "dT.ast") & dT.NF.melt$origin != "Other",], aes(origin, value, group = origin, colour = origin, fill = origin)) + geom_jitter(alpha = 0.5, size = 0.5) + coord_flip() + theme_classic() + scale_colour_viridis_d() + stat_compare_means(label.x = 0.5, label.y = 0.5) + theme(line = element_line(size = 1))
dT.origin.jitter <- dT.origin.jitter + theme(line = element_line(size = 1), axis.title.y = element_blank(), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10)) + ylab("dT")

# Density plot of ancestral dT values, grouped by strict 2R or 3R
# Jitter plot of ancestral dT values, grouped by strict 2R or 3R
# 2R vs 3R is NOT significantly different, Wilcoxon p-value = 0.14
redun.origin.density <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("dT.zeb.ancestral", "dT.ast.ancestral") & dT.NF.melt$origin != "Other",], aes(x = value, group = origin, fill = origin, colour = origin)) + geom_density(alpha = 0.1, size = 1) + theme_classic() + scale_colour_viridis_d() + scale_fill_viridis_d() + theme(legend.position = c(0.65,0.75), legend.title = element_blank())
redun.origin.density <- redun.origin.density + ylab("Density") + xlab("log(Ancestral dT)") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), line = element_line(size = 1)) + ylab("")

redun.origin.jitter <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("dT.zeb.ancestral", "dT.ast.ancestral") & dT.NF.melt$origin != "Other",], aes(origin, value, group = origin, colour = origin, fill = origin)) + geom_jitter(alpha = 0.5, size = 0.5) + coord_flip() + theme_classic() + scale_colour_viridis_d() + stat_compare_means(label.x = 0.5, label.y = 0.5) + theme(line = element_line(size = 1))
redun.origin.jitter <- redun.origin.jitter + theme(line = element_line(size = 1), axis.title.y = element_blank(), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10)) + ylab("Redundancy score") + xlab("")

# NEED TO CHANGE THE NF SCORE TO Specialization SCORE + CHANGE HOW IT IS CALCUATED (just # of new cell types) - then log(x+1) should work ~ok
# Density plot of specialization scores, grouped by strict 2R or 3R
# Jitter plot of specialization scores, grouped by strict 2R or 3R
# 2R vs 3R is NOT significantly different, Wilcoxon p-value = 0.22
sc.origin.density <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("NF.gene1.zeb", "NF.gene2.zeb", "NF.gene1.ast", "NF.gene2.ast") & dT.NF.melt$origin != "Other",], aes(x = log(value+1), group = origin, fill = origin, colour = origin)) + geom_density(alpha = 0.1, size = 1) + theme_classic() + scale_colour_viridis_d() + scale_fill_viridis_d() + theme(legend.position = c(0.65,0.75), legend.title = element_blank())
sc.origin.density <- sc.origin.density + ylab("Density") + xlab("log(Ancestral dT)") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank(), line = element_line(size = 1)) + ylab("")

sc.origin.jitter <- ggplot(dT.NF.melt[dT.NF.melt$variable == c("NF.gene1.zeb", "NF.gene2.zeb", "NF.gene1.ast", "NF.gene2.ast") & dT.NF.melt$origin != "Other",], aes(origin, log(value+1), group = origin, colour = origin, fill = origin)) + geom_jitter(alpha = 0.5, size = 0.5) + coord_flip() + theme_classic() + scale_colour_viridis_d() + stat_compare_means(label.x = 0.5, label.y = 0.5) 
sc.origin.jitter <- sc.origin.jitter + theme(line = element_line(size = 1), axis.title.y = element_blank(), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10)) + ylab("Specialization score") + xlab("")


# # Combine all 2R vs 3R and branch plots into 1 figure
# dev.new()
# dt.ecdf + redun.ecdf + sc.ecdf + dT.origin.density + redun.origin.density + sc.origin.density + dT.origin.jitter + redun.origin.jitter + sc.origin.jitter + plot_layout(nrow = 3, height = unit(c(40,30,10), c("mm", "mm")), width = unit(c(40), c("mm")), guides = "collect")


### Make DotPlots for divergent gene pairs

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

dot.zeb <- DotPlot(hypo.integrated.zeb, features = c("etv5b", "etv5a"), group.by = "integrated_Cluster", scale.max = 40) + RotatedAxis() + coord_flip() + scale_color_viridis() + theme(axis.title = element_blank(), axis.text = element_text(size = 8), axis.text.x = element_blank())
dot.ast <- DotPlot(hypo.integrated.ast, features = c("etv5b", "etv5a"), group.by = "integrated_Cluster", scale.max = 40) + RotatedAxis() + coord_flip() + scale_color_viridis() + theme(axis.title = element_blank(), axis.text = element_text(size = 8))

# Figure S4d
dots <- dot.zeb / dot.ast + plot_layout(nrow = 2, guides = "collect")



# Create matrices of normalized expression for each gene across integrated_SubclusterTypes

# Create lists of pairs divided by branch
library(plyr)
dTzeb.2 <- plyr::match_df(dT.list[[1]], dT.list[[2]], on = c("gene1", "gene2"))
dTast.2 <- plyr::match_df(dT.list[[2]], dT.list[[1]], on = c("gene1", "gene2"))
dT.combined <- merge(dTast.2, dTzeb.2, by = c("gene1", "gene2"))

dT.combined$genes <- paste(dT.combined$gene1, dT.combined$gene2, sep = "_")
dT.combined$origin <- ifelse(dT.combined$genes %in% twor$genes, "2R", ifelse(dT.combined$genes %in% threer$genes, "3R", "Other"))
dT.combined$branch <- dT.combined$branch.x

dT.combined$branch <- factor(dT.combined$branch, levels = branch.levels)
dT.combined$origin <- factor(dT.combined$origin, levels = c("Other", "3R", "2R"))

dTcombined <- lapply(unique(dT.combined$branch.x), function(x) dT.combined[dT.combined$branch.x == x,])
names(dTcombined) <- unique(dT.combined$branch.x)

# Load trinarization scores for cell types

trinarized.genes <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,".rds", sep = ""))

# Combine subcluster trinarization scores into dataframes by species

zeb.int <- Reduce(cbind, trinarized.genes$subcluster.zebrafish)
colnames(zeb.int) <- names(trinarized.genes$subcluster.zebrafish)
# row.names(zeb.int) <- row.names(trinarized.exp$subcluster.zebrafish[[1]])

ast.int <- Reduce(cbind, trinarized.genes$subcluster.astyanax)
colnames(ast.int) <- names(trinarized.genes$subcluster.astyanax)
# row.names(ast.int) <- row.names(trinarized.exp$subcluster.astyanax[[1]])

zeb.int.b <- ifelse(zeb.int > .95, 1, 0)
ast.int.b <- ifelse(ast.int > .95, 1, 0)

# # Match the gene pairs from the dTcombined lists from each .int matrix and cor() them
# 
# correlations <- lapply(dTcombined, function(x) apply(x, 1, function(x) cor(unlist(c(zeb.int[as.character(x[1]),intersect(colnames(zeb.int), colnames(ast.int))], zeb.int[as.character(x[2]),intersect(colnames(zeb.int), colnames(ast.int))])), unlist(c(ast.int[as.character(x[1]),intersect(colnames(zeb.int), colnames(ast.int))], ast.int[as.character(x[2]),intersect(colnames(zeb.int), colnames(ast.int))])))))
# correlations2 <- data.frame(branch = c(rep(names(correlations)[[1]], length(correlations[[1]])), rep(names(correlations)[[2]], length(correlations[[2]])), rep(names(correlations)[[3]], length(correlations[[3]])), rep(names(correlations)[[4]], length(correlations[[4]])), rep(names(correlations)[[5]], length(correlations[[5]])), rep(names(correlations)[[6]], length(correlations[[6]])), rep(names(correlations)[[7]], length(correlations[[7]])), rep(names(correlations)[[8]], length(correlations[[8]]))), value = unlist(correlations))
# correlations2$branch <- factor(correlations2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Danio rerio", "Astyanax mexicanus"))
# anno_df <- compare_means(value ~ branch, data = correlations2)
# 
# # Figure S7a
# correlation.trin.scores <- ggplot(correlations2, aes(x = branch, y = value)) + geom_jitter(color = "grey65", shape = 21) + geom_boxplot(fill = "transparent", outlier.color = NA) + theme_classic() + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + stat_compare_means()
# 
## Do it for binarized expression (there or not)

correlations.b <- sapply(dTcombined, function(x) apply(x, 1, function(x) cor(unlist(c(zeb.int.b[as.character(x[1]),intersect(colnames(zeb.int.b), colnames(ast.int.b))], zeb.int.b[as.character(x[2]),intersect(colnames(zeb.int.b), colnames(ast.int.b))])), unlist(c(ast.int.b[as.character(x[1]),intersect(colnames(zeb.int.b), colnames(ast.int.b))], ast.int.b[as.character(x[2]),intersect(colnames(zeb.int.b), colnames(ast.int.b))])))))

correlations.b <- lapply(seq_along(correlations.b), function(x) {
  names(correlations.b[[x]]) <- dTcombined[[x]]$genes
  return(correlations.b[[x]])
  })


correlations.b.2 <- data.frame(branch = c(rep(names(correlations)[[1]], length(correlations.b[[1]])), rep(names(correlations)[[2]], length(correlations.b[[2]])), rep(names(correlations)[[3]], length(correlations.b[[3]])), rep(names(correlations)[[4]], length(correlations.b[[4]])), rep(names(correlations)[[5]], length(correlations.b[[5]])), rep(names(correlations)[[6]], length(correlations.b[[6]])), rep(names(correlations)[[7]], length(correlations.b[[7]])), rep(names(correlations)[[8]], length(correlations.b[[8]]))), value = unlist(correlations.b))
correlations.b.2$branch <- factor(correlations.b.2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Danio rerio", "Astyanax mexicanus"))
correlations.b.2$origin <- ifelse(rownames(correlations.b.2) %in% twor$genes, "2R", ifelse(dT.combined$genes %in% threer$genes, "3R", "Other"))


anno_df <- compare_means(value ~ branch, data = correlations.b.2)

# Figure S7b
correlation.exp <- ggplot(correlations.b.2[correlations.b.2$branch != "Opisthokonta",], aes(x = value, y = branch, colour = branch, fill = branch, group = branch)) + geom_density_ridges(alpha = 0.1, scale = 2, jittered_points = T, point_size = 0.5) 
correlation.exp <- correlation.exp + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme_classic() + ylab("") + xlab("")
correlation.origin <- ggplot(correlations.b.2[correlations.b.2$branch != "Opisthokonta" & correlations.b.2$origin != "Other",], aes(x = value, y = origin, colour = origin, fill = origin)) + geom_density_ridges(alpha = 0.25, scale = 2, jittered_points = T, point_size = 0.5) 
correlation.origin <- correlation.origin + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme_classic() + scale_colour_viridis_d() + scale_fill_viridis_d() + ylab("") + xlab("Paralog expression correlation")


# Combine all 2R vs 3R and branch plots into 1 figure
dev.new()
dt.ecdf + redun.ecdf + sc.ecdf + dT.origin.density + redun.origin.density + sc.origin.density + dT.origin.jitter + redun.origin.jitter + sc.origin.jitter + plot_layout(ncol = 3, height = unit(c(40,30,10), c("mm", "mm")), width = unit(c(47), c("mm")), guides = "collect")

design2 <- "
ABD
CBD"

dev.new()
dot.zeb + correlation.exp + dot.ast + correlation.origin + plot_layout(nrow = 10, ncol = 5, height = unit(c(10), c("mm")), width = unit(c(80,30,30), c("mm")), design = design2, guides = "collect")




