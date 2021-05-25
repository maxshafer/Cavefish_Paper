library(reshape2)
library(Seurat)
library(viridis)
library(ggrepel)
library(scales)
library(ggsignif)
library(ggpubr)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")

# Load GENIE3 files from all the comps, for cave and surface

gene.lists <- c("nps", "nts", "synaptic", "ion")

surface.lists <- lapply(gene.lists, function(x) readRDS(paste("surface/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
cave.lists <- lapply(gene.lists, function(x) readRDS(paste("cave/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
surface.modules <- lapply(gene.lists, function(x) readRDS(paste("surface/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
cave.modules <- lapply(gene.lists, function(x) readRDS(paste("cave/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))

names(surface.lists) <- gene.lists
names(cave.lists) <- gene.lists
names(surface.modules) <- gene.lists
names(cave.modules) <- gene.lists

## Make functions for extracting/plotting data
## load functions
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/SCENIC_functions_new.R")

## Load GO_lists and make plot data and plots

GO_lists <- lapply(gene.lists, function(x) union(unique(surface.lists[[x]]$Target), unique(cave.lists[[x]]$Target)))
names(GO_lists) <- gene.lists
GO_data <- lapply(seq_along(GO_lists), function(y) lapply(seq_along(GO_lists[[y]]), function(x) extractPlotData(surface = surface.lists[[y]], cave = cave.lists[[y]], cor.surface = surface.modules[[y]], cor.cave = cave.modules[[y]], np = GO_lists[[y]][x])))

names(GO_data) <- names(GO_lists)
for (i in 1:length(GO_data)) {
	names(GO_data[[i]]) <- GO_lists[[i]]
}

GO_plots <- lapply(seq_along(GO_lists), function(y) lapply(names(GO_data[[y]]), function(x) plotPlotData(data = GO_data[[y]], gene = x, quantile = 0.98, log = F)))

names(GO_plots) <- names(GO_lists)
for (i in 1:length(GO_plots)) {
	names(GO_plots[[i]]) <- GO_lists[[i]]
}

# GO plots

galn <- GO_plots[["nps"]][["galn"]] + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10, colour = "black"), legend.position = "none")
hcrt <- GO_plots[["nps"]][["hcrt"]] + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10, colour = "black"), legend.position = "none")
oxt <- GO_plots[["nps"]][["oxt"]] + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10, colour = "black"), legend.position = "none")
avp <- GO_plots[["nps"]][["avp"]] + theme(axis.text = element_text(size = 8, colour = "black"), axis.title = element_text(size = 10, colour = "black"), legend.position = "none")

# GO_plots[["nps"]][["hcrt"]]


## Calculate the Similarity Index for each NP for plotting
## Use weight cutoffs for calling GRNs, 0.001 and 0.005
## So tfModules_asDF.surface and tfModules_asDF.cave have done this for me already

surface.grns <- lapply(surface.modules, function(x) x[x$method == "top50",])

for (i in 1:length(surface.grns)) {
	names <- unique(surface.grns[[i]]$Target)
	surface.grns[[i]] <- lapply(names, function(x) surface.grns[[i]][surface.grns[[i]]$Target == x, "TF"])
	names(surface.grns[[i]]) <- names
}

cave.grns <- lapply(cave.modules, function(x) x[x$method == "top50",])

for (i in 1:length(cave.grns)) {
	names <- unique(cave.grns[[i]]$Target)
	cave.grns[[i]] <- lapply(names, function(x) cave.grns[[i]][cave.grns[[i]]$Target == x, "TF"])
	names(cave.grns[[i]]) <- names
}


for (i in 1:length(surface.grns)) {
	surface.grns[[i]] <- surface.grns[[i]][intersect(names(surface.grns[[i]]), names(cave.grns[[i]]))]
	cave.grns[[i]] <- cave.grns[[i]][intersect(names(surface.grns[[i]]), names(cave.grns[[i]]))]
}

conserved.grns <- lapply(seq_along(cave.grns), function(y) lapply(seq_along(cave.grns[[y]]), function(x) intersect(surface.grns[[y]][[x]], cave.grns[[y]][[x]])))

names(conserved.grns) <- names(cave.grns)
for (i in 1:length(conserved.grns)) {
	names(conserved.grns[[i]]) <- names(cave.grns[[i]])
}


## Calculate the Similarity Index

SI <- lapply(seq_along(conserved.grns), function(x) calcSimilarityIndex(conserved = conserved.grns[[x]], species.1 = surface.grns[[x]], species.2 = cave.grns[[x]]))

names(SI) <- names(cave.grns)
for (i in 1:length(SI)) {
	names(SI[[i]]) <- names(cave.grns[[i]])
}

## Calculate the % paralogs in each GRN
### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("../../Seurat_v3_Integration/mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("../../Seurat_v3_Integration/mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())]#[3:8]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))

names(go_lists) <- list.files()[grep("GO", list.files())]#[3:8]
names(go_lists)


paralog.numbers <- list()
paralog.numbers <- lapply(seq_along(conserved.grns), function(y) lapply(seq_along(conserved.grns[[y]]), function(x) calcParalog(conserved = conserved.grns[[y]], species.1 = surface.grns[[y]], species.2 = cave.grns[[y]], ngenes.1 = length(go_lists[[13]][go_lists[[13]] %in% row.names(GetAssayData(hypo.ast))]), ngenes.2 = length(go_lists[[13]][go_lists[[13]] %in% row.names(GetAssayData(hypo.ast))]), i = x)))

grn.paralog <- lapply(paralog.numbers, function(x) data.frame(do.call(rbind, x)))
grn.paralog <- lapply(grn.paralog, function(x) rbind(x, colSums(x)))

for (i in 1:length(grn.paralog)) {
	row.names(grn.paralog[[i]]) <- c(names(conserved.grns[[i]]), "sums")
	grn.paralog[[i]]$percent.para <- unlist(apply(grn.paralog[[i]], 1, function(x) sum(x[1], x[5])/sum(x[1], x[5], x[3], x[7]))*100)
	grn.paralog[[i]]$percent.conserved.1 <- unlist(apply(grn.paralog[[i]], 1, function(x) x[9]/sum(x[9], x[10])*100))
	grn.paralog[[i]]$percent.conserved.2 <- unlist(apply(grn.paralog[[i]], 1, function(x) x[9]/sum(x[9], x[11])*100))
	# grn.paralog[[i]]$fisher.pval <- as.numeric(fisher.results[,1])
	grn.paralog[[i]]$fisher.odds.1 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[1:4] + x[5:8], ncol = 2, nrow = 2))))[,3])
	grn.paralog[[i]]$fisher.odds.pval.1 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[1:4] + x[5:8], ncol = 2, nrow = 2))))[,1])
	#grn.paralog[[i]]$fisher.odds.2 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2))))[,3])
	#grn.paralog[[i]]$fisher.odds.pval.2 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2))))[,1])
	grn.paralog[[i]]$fisher.odds.3 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2))))[,3])
	grn.paralog[[i]]$fisher.odds.pval.3 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2))))[,1])
	grn.paralog[[i]]$gene <- row.names(grn.paralog[[i]])
}


## Calculate the correlation between the GENIE3 results for each species, for each GOI

correlations <- lapply(seq_along(GO_lists), function(y) unlist(lapply(seq_along(GO_lists[[y]]), function(x) cor(GO_data[[y]][[x]]$weight_cave, GO_data[[y]][[x]]$weight_surface))))

names(correlations) <- names(GO_lists)
for (i in 1:length(correlations)) {
	names(correlations[[i]]) <- GO_lists[[i]]
}

## Put Correlation and SI into the same df
## gene / GO / SI / corr

SI <- lapply(SI, function(x) unlist(x))
correlations <- lapply(correlations, function(x) unlist(x))
correlations <- lapply(seq_along(correlations), function(x) correlations[[x]][names(SI[[x]])])
names(correlations) <- names(SI)

df <- reshape2::melt(do.call(rbind, lapply(seq_along(SI), function(x) data.frame(gene = names(SI[[x]]), GO = names(SI)[x], SI = SI[[x]], corr = correlations[[x]], percent.para = grn.paralog[[x]]$percent.para[1:length(SI[[x]])]))))

df2 <- do.call(rbind, lapply(seq_along(SI), function(x) data.frame(gene = names(SI[[x]]), GO = names(SI)[x], fisher.odds = grn.paralog[[x]]$fisher.odds.1[1:length(SI[[x]])], fisher.odds.pval = grn.paralog[[x]]$fisher.odds.pval.1[1:length(SI[[x]])])))
df2$GO <- factor(df2$GO, levels = c("nps", "nts", "synaptic", "ion"))
# Re-order factor levels
# Some genes are in multiple categories, and therefore can't be factored correctly...

df <- df[!duplicated(df[c(1,3)]),]

df$gene <- factor(df$gene, levels = df[order(df$value[df$variable == "corr"]), "gene"], ordered = T)
df$variable <- factor(df$variable, levels = c("corr", "SI", "percent.para"))
df$GO <- factor(df$GO, levels = c("nps", "nts", "synaptic", "ion"))
## Plot scatter plots and geom_tile heatmaps!

anno_df <- compare_means(value ~ GO, group.by = 'variable', data = df)
anno_df$y_pos <- c( 0.75, 0.9, 0.95, 1.0, 1.05, 1.1, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 15, 16, 17, 18, 19, 20)

scatter1 <- ggplot(df[df$variable == "corr",], aes(x = GO, y = value, color = variable)) + geom_point(position = position_jitterdodge(), size = .5) + geom_boxplot(fill = NA, size = .5, outlier.shape = NA, color = "black") + scale_color_manual(values = c("grey55")) + theme_classic() + theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8)) + scale_x_discrete(labels = c("Neuropeptides", "Neurotransmitters", "Synaptic Genes", "Ion Channel")) + guides(color = F) #+ geom_signif(data = anno_df[anno_df$variable == "corr",], aes(xmin = group1, xmax = group2, annotations = p.signif, y_position = y_pos), manual = T, textsize = 3) + ylab("Correlation")
scatter2 <- ggplot(df[df$variable == "SI",], aes(x = GO, y = value, color = variable)) + geom_point(position = position_jitterdodge(), size = .5) + geom_boxplot(fill = NA, size = .5, outlier.shape = NA, color = "black") + scale_color_manual(values = c("red")) + theme_classic() + theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8)) + scale_x_discrete(labels = c("Neuropeptides", "Neurotransmitters", "Synaptic Genes", "Ion Channel")) + guides(color = F) #+ geom_signif(data = anno_df[anno_df$variable == "SI",], aes(xmin = group1, xmax = group2, annotations = p.signif, y_position = y_pos), manual = T, textsize = 3) + ylab("Similarity Index (SI)")
scatter3 <- ggplot(df[df$variable == "percent.para",], aes(x = GO, y = value, color = variable)) + geom_point(position = position_jitterdodge(), size = .5) + geom_boxplot(fill = NA, size = .5, outlier.shape = NA, color = "black") + scale_color_manual(values = c("black")) + theme_classic() + theme(axis.title.x = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8)) + scale_x_discrete(labels = c("Neuropeptides", "Neurotransmitters", "Synaptic Genes", "Ion Channel")) + guides(color = F)# + geom_signif(data = anno_df[anno_df$variable == "percent.para",], aes(xmin = group1, xmax = group2, annotations = p.signif, y_position = y_pos), manual = T, textsize = 3) + ylab("% Paralogs")

odds.scatter <- ggplot(df2, aes(x = GO, y = fisher.odds, color = ifelse(fisher.odds.pval < 0.05, "yes", "no"))) + geom_jitter(size = 0.5) + geom_boxplot(fill = NA, size = 0.5, outlier.shape = NA, color = "black") + scale_color_manual(values = c("black", "red")) + theme_classic() + theme(axis.title = element_text(size = 10), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8)) + scale_x_discrete(labels = c("Neuropeptides", "Neurotransmitters", "Synaptic Genes", "Ion Channel")) + ylab("Odds ratio") + xlab("") + guides(color = F)

scatter2 <- scatter2 + ylab("Similarity Index (SI)") + theme(text = element_text(colour = "black"))

## Make small figure for galnain and oxytocin cluster expression

galn.ob <- subset(hypo.ast, idents = "Neuronal_07")
oxt.ob <- subset(hypo.ast, idents = "Neuronal_19")

dot.galn <- DotPlot(galn.ob, features = c("galn"), group.by = "morph", scale.min = 30, scale.max = 80, dot.scale = 4) 
dot.galn <- dot.galn + theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8)) + scale_color_viridis_c(limits = c(-1.5,1.5), option = "B")

dot.oxt <- DotPlot(oxt.ob, features = c("oxt", "avp", "ENSAMXG00000021172"), group.by = "morph", scale.min = 30, scale.max = 80, dot.scale = 4) 
dot.oxt <- dot.oxt + theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank()) + scale_color_viridis_c(limits = c(-1.5,1.5), option = "B")

dots <- dot.galn + dot.oxt + plot_layout(guides = "collect", widths = unit(c(7,13), "mm"), height = unit(c(35,35), "mm"))

# Put them together

design <- "
ABDE
CCFG"

pdf("../Figures/Hypo_AstMex_TF-analysis-figures.pdf", height = 11, width = 9) 
dots + scatter2 + galn + hcrt + oxt + avp + plot_layout(design = design, height = unit(c(30), c("mm")), width = unit(c(10,20,30,30), c("mm")))
dev.off()



