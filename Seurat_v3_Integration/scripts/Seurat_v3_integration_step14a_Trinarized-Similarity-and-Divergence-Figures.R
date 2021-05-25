library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(viridis)
library(ggpubr)
library(venneuler)


# This script makes figures which show that trinarized genes (all genes expressed by each subcluster) are highly similar both for orthologous, but also 'paralogous' subclusters
# This is due to a large signal from housekeeping genes, which are also enriched for paralogs like subcluster specific genes

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

a = 1.5
b = 2
f = 0.1

## Load trinarized gene lists

# trinarized.genes <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,".rds", sep = ""))
trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# Load hk filtered file
trinarized.exp.hkfilt <- readRDS(file = paste("trinarized_expression_hk-filtered_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# Add species-specific subclusters
norm.cluster.filtered.zeb <- lapply(c(trinarized.exp$subcluster.zebrafish, trinarized.exp$specific.zebrafish), function(x) names(x))
norm.cluster.filtered.ast <- lapply(c(trinarized.exp$subcluster.astyanax, trinarized.exp$specific.astyanax), function(x) names(x))


SI.list <- readRDS(paste("SI_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))
SI.list.hkfilt <- readRDS(paste("SI_trinarized_hk-filtered_a",a,"_b",b, "_f",f,".rds", sep = ""))

# SI <- SI.list[[1]]
# SI.sub <- SI.list[[2]]
# SI.ast <- SI.list[[3]]
# SI.ast.sub <- SI.list[[4]]
# 
# colnames(SI.sub) <- c("Subcluster", "values", "Cluster")

## Similarity Index (SI) and % paralog figures

matrix.ast2.zeb.sub <- SI.list[[11]]
names(matrix.ast2.zeb.sub) <- names(matrix.ast2.zeb.sub[[1]])


## Make row-scaled heatmap
matrix.data <- lapply(matrix.ast2.zeb.sub, function(x) unlist(x))
matrix.data <- reshape2::melt(lapply(matrix.data, function(x) (x)/max(x)))
# matrix.data <- reshape2::melt(matrix.data)
matrix.data$L2 <- rep(names(matrix.ast2.zeb.sub[[1]]),151)
matrix.data$L1 <- factor(matrix.data$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
matrix.data$L2 <- factor(matrix.data$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
SI.matrix.scaled <- ggplot(matrix.data, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6)) + scale_fill_viridis(direction = 1, limits = c(0,1))

SI.matrix.scaled <- SI.matrix.scaled + theme(axis.title = element_blank())


## Make comparison figure SI between zeb-ast and cave-surface before and after HK filtering

SI.combined.sub <- rbind(SI.list[[2]], SI.list[[4]])
SI.combined.sub.hkfilt <- rbind(SI.list.hkfilt[[2]], SI.list.hkfilt[[4]])
SI.sub <- SI.combined.sub
SI.sub$hkfilt <- SI.combined.sub.hkfilt$values
SI.sub$species <- c(rep("Between species", 151), rep("Between species-morphs", 151))
colnames(SI.sub) <- c("Subcluster", "trinarized_SI", "Cluster", "trinarized_SI_hk-filtered", "Comparison")
SI.sub <- reshape2::melt(SI.sub)

SI.plot <- ggplot(SI.sub, aes(x = interaction(variable, Comparison), y = value, colour = Comparison, group = interaction(Comparison, variable))) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic() + scale_colour_viridis_d()
SI.plot <- SI.plot + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")"))) + NoLegend()

## Make Correlation matrix

norm.cluster <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Normed_expression_data.rds")
str(norm.cluster, max.level = 2)

intersect.genes <- intersect(row.names(norm.cluster[[4]][[1]]), row.names(norm.cluster[[4]][[2]]))

zeb.avg <- norm.cluster[[4]][[1]][intersect.genes,names(matrix.ast2.zeb.sub[[1]])]
ast.avg <- norm.cluster[[4]][[2]][intersect.genes,names(matrix.ast2.zeb.sub[[1]])]

zeb.ast <- cor(zeb.avg, ast.avg, method = "pearson") # row x column
# zeb.ast <- zeb.ast/rowMax(zeb.ast) # scale by row, looks very similar
zeb.ast <- reshape2::melt(zeb.ast)
zeb.ast$Var1 <- factor(zeb.ast$Var1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
zeb.ast$Var2 <- factor(zeb.ast$Var2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))

Corr.matrix.scaled <- ggplot(zeb.ast, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 6)) + scale_fill_viridis(direction = 1, limits = c(0,1))
Corr.matrix.scaled <- Corr.matrix.scaled + theme(axis.title = element_blank())


## Make comparison of SI for between species and species-morphs

SI.list <- readRDS(paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

SI.combined.sub <- rbind(SI.list[[2]], SI.list[[4]])
SI.combined.sub$species <- c(rep("Between species", nrow(SI.list[[2]])), rep("Between species-morphs", nrow(SI.list[[4]])))
colnames(SI.combined.sub) <- c("Subcluster", "SI", "Cluster", "Comparison")

species.SI.sub <- ggplot(SI.combined.sub, aes(x = Comparison, y = SI, colour = Comparison, group = Comparison)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic() + scale_colour_viridis_d()
species.SI.sub <- species.SI.sub + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")"))) + NoLegend()


## Combine the above to make first Supplemental Figure
## Put them together (minus the venn diagram, which get's a plot spacer)

design <- "
AAAABBBB
CCCCEEDD"

pdf("Figures/Hypo_integrated_corr_trinarization_hk_plots.pdf", height = 11, width = 9) 
(Corr.matrix.scaled + SI.matrix.scaled + plot_layout(ncol = 2, widths = unit(c(60), c("mm")), height = unit(c(60), c("mm")), guides = "collect")) / (SI.plot + plot_spacer() + species.SI.sub + plot_layout(nrow = 1, height = unit(c(30), c("mm")), width = unit(c(45,60,15), c("mm")), guides = "collect")) + plot_layout(nrow = 2)
dev.off()


## Then plot the venn diagram, and combine in Affinity/Illustrator

hk <- readRDS(file = paste("housekeeping_genes_a",a,"_b",b, "_f",f,".rds", sep = ""))

## Calculate housekeeping SI

Gt <- intersect(names(hk[[1]][hk[[1]]]), names(hk[[2]][hk[[2]]]))
Ga <- names(hk[[1]][hk[[1]]])
Gb <- names(hk[[2]][hk[[2]]])

hk.SI <- 1 - sqrt( ( 1 - length(Gt)/length(Ga)) * (1 - length(Gt)/length(Gb)) )


## Calculate the % of paralogs for HK gene sets and make Venn diagram

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)

mart.1 <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart.2 <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)

conserved <- Gt
genes.1 <- Ga[!(Ga %in% Gt)]
genes.2 <- Gb[!(Gb %in% Gt)]

paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(conserved, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(conserved, mart.2$Gene.name)])
paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart.1$Gene.name)])
paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart.2$Gene.name)])
paralog.union <- union(paralog.con, union(paralog.1, paralog.2))

paralog.1 <- paralog.1[paralog.1 %in% mart.2$Gene.name]
paralog.2 <- paralog.2[paralog.2 %in% mart.1$Gene.name]

a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.2)])
b1 <- length(genes.1) - a1
c1 <- length(union(paralog.con, paralog.2)) - a1
d1 <- length(genes.1)- b1

a2 <- length(genes.2[genes.2 %in% union(paralog.con, paralog.1)])
b2 <- length(genes.2) - a2
c2 <- length(union(paralog.con, paralog.1)) - a2
d2 <- length(genes.2) - b2

vec <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2)

percent.para.1 <- vec[1]/sum(vec[1], vec[3])*100
percent.para.2 <- vec[5]/sum(vec[5], vec[7])*100

fisher.results.1 <- fisher.test(matrix(data = vec[1:4], ncol = 2, nrow = 2)) # Zebrafish paralog enrichment
fisher.results.2 <- fisher.test(matrix(data = vec[5:8], ncol = 2, nrow = 2)) # Mexican tetra paralog enrichment


## Make a venn diagram which shows the overlap, as well as the % paralogs (as sub circles of the species-specific HK genes)

venn.data <- c(Zebrafish_HK = length(genes.1), 
               Mexican_tetra_HK = length(genes.2),
               Zebrafish_HK_paralogs = 0,
               Mexican_tetra_HK_paralogs = 0,
               'Zebrafish_HK&Mexican_tetra_HK' = length(Gt),
               'Zebrafish_HK&Zebrafish_HK_paralogs' = as.numeric(vec[1]),
               'Zebrafish_HK&Mexican_tetra_HK_paralogs' = 0,
               'Mexican_tetra_HK&Zebrafish_HK_paralogs' = 0,
               'Mexican_tetra_HK&Mexican_tetra_HK_paralogs' = as.numeric(vec[5]),
               'Zebrafish_HK_paralogs&Mexican_tetra_HK_paralogs' = 0)

vd1 <- venneuler(venn.data)

# For some reason, this doesn't make the HK paralog circles completely overlapping with the HK circles
dev.new()
plot(vd1)















# ## Plot progenitors versus neurons
# 
# SI.sub$Neuronal <- "Neuronal"
# SI.sub$Neuronal[grepl("Progenitors", SI.sub$Subcluster)] <- "Progenitors"
# SI.sub$Neuronal <- factor(SI.sub$Neuronal, levels = c("Progenitors", "Neuronal"))
# neuronal <- ggplot(as.data.frame(SI.sub), aes(x = Neuronal, y = values, group = Neuronal, color = Neuronal)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic()
# neuronal <- neuronal + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")")))
# 
# 
# SI.compare <- data.frame(cell_type = SI$Cluster, cluster = SI$values, cell_types = aggregate(SI.sub[,2], list(SI.sub$Cluster), mean)[,2])
# cluster <- ggplot(reshape::melt(SI.compare), aes(x = variable, y = value, group = cell_type, color = cell_type)) + geom_jitter(size = 0.5) + geom_line(aes(group = cell_type), colour = "grey75") + geom_boxplot(aes(group = variable),outlier.color = NA, fill = "transparent", color = "black") + theme_classic() + theme(legend.position = "none")
# cluster <- cluster + scale_x_discrete(labels=c("cluster" = "Cluster", "cell_types" = "Subcluster")) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab(expression(paste("Similarity Index (", italic("SI"), ")"))) 
# 
# t.test(SI.sub[SI.sub$Neuronal == "Neuronal", "values"], SI.sub[SI.sub$Neuronal == "Progenitors", "values"])
# t.test(SI.compare$cell_types, SI.compare$cluster, paired = T)
# 
# 
# ### Paralog plots
# ## Make SI by Paralog (sum?? or mean??)
# 
# para.results <- readRDS(paste("Paralog-results-trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))
# # para.results <- readRDS(paste("Paralog-results-trinarized_hk-filtered_a",a,"_b",b, "_f",f,".rds", sep = ""))
# 
# para.data <- reshape2::melt(para.results[[1]][,13:18])
# para.data.sub <- reshape2::melt(para.results[[2]][,13:19])
# 
# percent.para <- ggplot(para.data.sub[para.data.sub$variable == c("percent.para.1", "percent.para.2"),], aes(x = variable, y = value, colour = variable, group = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic() + ylim(c(0,18))
# percent.para <- percent.para + scale_x_discrete(labels=c("percent.para.1" = "Zebrafish", "percent.para.2" = "Mexican tetra")) + theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) + ylab("% Paralogs")
# 
# odds.para <- ggplot(para.data.sub[para.data.sub$variable == c("fisher.odds.1", "fisher.odds.2"),], aes(x = variable, y = value, colour = variable, group = variable)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + theme_classic()
# odds.para <- odds.para + scale_x_discrete(labels=c("fisher.odds.1" = "Zebrafish", "fisher.odds.2" = "Mexican tetra")) + theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_blank()) + ylab("Odds Ratio")
# 
# # Put together Figure 2
# 
# part1 <- SI.matrix.scaled + (percent.para + odds.para + plot_layout(ncol = 2)) + plot_layout(nrow = 2, heights = c(3,1), guides = "collect")
# part1 <- part1 + plot_layout(height = unit(c(90,30), c("mm", "mm")), width = unit(c(90), c("mm")))
# dev.new()
# part1
# 
# 
# 
# 
# ## Make paralog % by drift scatter
# 
# para.results[[2]]$Neuronal <- "Non-Neuronal"
# para.results[[2]]$Neuronal[grepl("Neuronal", para.results[[2]]$Subcluster)] <- "Neuronal"
# 
# r <- cor((para.results[[2]]$percent.para.2 + para.results[[2]]$percent.para.1)/2, para.results[[2]]$similarityindex)
# para.results.neuronal <- para.results[[2]][para.results[[2]]$Neuronal == "Neuronal",]
# r <- cor((para.results.neuronal$percent.para.2 + para.results.neuronal$percent.para.1)/2, para.results.neuronal$similarityindex)
# para.results.non <- para.results[[2]][para.results[[2]]$Neuronal == "Non-Neuronal",]
# r <- cor((para.results.non$percent.para.2 + para.results.non$percent.para.1)/2, para.results.non$similarityindex)
# 
# ## If you group by neuronal/non, there is a clear difference in the relationship
# para.drift <- ggplot(para.results[[2]], aes(x = (percent.para.2 + percent.para.1)/2, y = similarityindex, color = Cluster, label = Subclusters)) + geom_point(size = 0.5) + guides(color = F) + ylab("Drift Index") + xlab("Mean of % Paralogs per species") # + geom_smooth(method = lm, se = T, color = "black", fill = "grey85")
# para.drift <- para.drift + theme_classic() + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + annotate(geom = "text", x = 2.5, y = 0.35, label = paste("Pearson r = ", round(r,digits = 2), sep = ""), color = "black") + NoLegend()
# 
# 
# ## Row-scaled delta between SI and Corrected-SI
# SI.delta <- lapply(seq_along(SI.list[[11]]), function(x) lapply(seq_along(SI.list[[11]][[x]]), function(y) SI.list[[12]][[x]][[y]] - SI.list[[11]][[x]][[y]]))
# names(SI.delta) <- names(gene.lists.pos$subcluster.conserved)
# 
# SI.delta <- lapply(SI.delta, function(x) unlist(x))
# SI.delta <- reshape2::melt(lapply(SI.delta, function(x) (x)/max(x)))
# SI.delta$L2 <- rep(names(SI.list[[12]][[1]]),151)
# SI.delta$L1 <- factor(SI.delta$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# SI.delta$L2 <- factor(SI.delta$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# 
# SI.delta.plot <- ggplot(SI.delta, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 1.5), axis.title = element_blank()) + scale_fill_viridis(direction = 1, option = "C")
# 
# matrix.data <- lapply(SI.list[[11]], function(x) unlist(x))
# names(matrix.data) <- names(gene.lists.pos$subcluster.conserved)
# matrix.data <- reshape2::melt(matrix.data)
# matrix.data$L2 <- rep(names(matrix.ast2.zeb.sub[[1]]),151)
# matrix.data$L1 <- factor(matrix.data$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# matrix.data$L2 <- factor(matrix.data$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# 
# matrix.data.corr <- lapply(SI.list[[12]], function(x) unlist(x))
# names(matrix.data.corr) <- names(gene.lists.pos$subcluster.conserved)
# matrix.data.corr <- reshape2::melt(matrix.data.corr)
# matrix.data.corr$L2 <- rep(names(matrix.ast2.zeb.sub[[1]]),151)
# matrix.data.corr$L1 <- factor(matrix.data.corr$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# matrix.data.corr$L2 <- factor(matrix.data.corr$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# 
# SI.combined.corr <- rbind(matrix.data[matrix.data$L1 == matrix.data$L2,], matrix.data.corr[matrix.data.corr$L1 == matrix.data.corr$L2,])
# SI.combined.corr$sample <- c(rep("Uncorrected SI", 151), rep("Corrected SI", 151))
# SI.combined.corr$sample <- factor(SI.combined.corr$sample, levels = c("Uncorrected SI", "Corrected SI"))
# 
# corr.delta.plot <- ggplot(SI.combined.corr, aes(x = sample, y = value, color = sample)) + geom_jitter(size = 0.5) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d(option = "C") + scale_color_manual(values = c("#5D01A6FF", "#FDB32FFF")) + theme_classic()
# corr.delta.plot <- corr.delta.plot + theme(axis.text = element_text(size = 8), axis.title.x = element_blank())+ ylab("Δ SI vs Corrected SI") + NoLegend()
# 
# # Put them together into S5
# 
# # dev.new()
# # corr.delta.plot + plot_layout(nrow = 1, widths = unit(c(10), "mm"), height = unit(c(20), "mm"), guides = "collect")
# # dev.new()
# # para.drift + plot_layout(nrow = 1, widths = unit(c(25), "mm"), height = unit(c(25), "mm"), guides = "collect")
# # dev.new()
# # SI.delta.plot + plot_layout(widths = unit(c(90), "mm"), height = unit(c(90), "mm"), guides = "collect")
# 
# part2 <- SI.delta.plot + (para.drift + corr.delta.plot + plot_layout(ncol = 2)) + plot_layout(nrow = 2, heights = c(3,1), guides = "collect")
# part2 <- part2 + plot_layout(height = unit(c(90,30), c("mm", "mm")), width = unit(c(90), c("mm")), guides = "collect")
# dev.new()
# part2
# 
# 
# 
# 
# 
