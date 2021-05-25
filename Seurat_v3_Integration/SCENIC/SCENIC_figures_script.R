library(reshape2)
library(Seurat)
library(viridis)
library(ggrepel)
library(scales)
library(ggsignif)
library(ggpubr)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC")
hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vR.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")

# Load GENIE3 files from all the comps, for ast and zeb

gene.lists <- c("nps", "nts", "synaptic", "ion")

zeb.lists <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
ast.lists <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
zeb.modules <- lapply(gene.lists, function(x) readRDS(paste("drerio/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
ast.modules <- lapply(gene.lists, function(x) readRDS(paste("amexicanus/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))

names(zeb.lists) <- gene.lists
names(ast.lists) <- gene.lists
names(zeb.modules) <- gene.lists
names(ast.modules) <- gene.lists

## Make functions for extracting/plotting data
## load functions
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/SCENIC_functions.R")

## Load GO_lists and make plot data and plots

GO_lists <- lapply(gene.lists, function(x) union(unique(zeb.lists[[x]]$Target), unique(ast.lists[[x]]$Target)))
names(GO_lists) <- gene.lists
GO_data <- lapply(seq_along(GO_lists), function(y) lapply(seq_along(GO_lists[[y]]), function(x) extractPlotData(zeb = zeb.lists[[y]], ast = ast.lists[[y]], cor.zeb = zeb.modules[[y]], cor.ast = ast.modules[[y]], np = GO_lists[[y]][x])))

names(GO_data) <- names(GO_lists)
for (i in 1:length(GO_data)) {
	names(GO_data[[i]]) <- GO_lists[[i]]
}

GO_plots <- lapply(seq_along(GO_lists), function(y) lapply(names(GO_data[[y]]), function(x) plotPlotData(data = GO_data[[y]], gene = x, quantile = 0.98, log = F)))

names(GO_plots) <- names(GO_lists)
for (i in 1:length(GO_plots)) {
	names(GO_plots[[i]]) <- GO_lists[[i]]
}

# GO_plots[["nps"]][["oxt"]]


## Calculate the Drift Index for each NP for plotting
## Use weight cutoffs for calling GRNs, 0.001 and 0.005
## So tfModules_asDF.zeb and tfModules_asDF.ast have done this for me already

zeb.grns <- lapply(zeb.modules, function(x) x[x$method == "top50",])

for (i in 1:length(zeb.grns)) {
	names <- unique(zeb.grns[[i]]$Target)
	zeb.grns[[i]] <- lapply(names, function(x) zeb.grns[[i]][zeb.grns[[i]]$Target == x, "TF"])
	names(zeb.grns[[i]]) <- names
}

ast.grns <- lapply(ast.modules, function(x) x[x$method == "top50",])

for (i in 1:length(ast.grns)) {
	names <- unique(ast.grns[[i]]$Target)
	ast.grns[[i]] <- lapply(names, function(x) ast.grns[[i]][ast.grns[[i]]$Target == x, "TF"])
	names(ast.grns[[i]]) <- names
}


for (i in 1:length(zeb.grns)) {
	zeb.grns[[i]] <- zeb.grns[[i]][intersect(names(zeb.grns[[i]]), names(ast.grns[[i]]))]
	ast.grns[[i]] <- ast.grns[[i]][intersect(names(zeb.grns[[i]]), names(ast.grns[[i]]))]
}

conserved.grns <- lapply(seq_along(ast.grns), function(y) lapply(seq_along(ast.grns[[y]]), function(x) intersect(zeb.grns[[y]][[x]], ast.grns[[y]][[x]])))

names(conserved.grns) <- names(ast.grns)
for (i in 1:length(conserved.grns)) {
	names(conserved.grns[[i]]) <- names(ast.grns[[i]])
}


## Calculate the Drift Index

SI <- lapply(seq_along(conserved.grns), function(x) calcDriftIndex(conserved = conserved.grns[[x]], species.1 = zeb.grns[[x]], species.2 = ast.grns[[x]]))

names(SI) <- names(ast.grns)
for (i in 1:length(SI)) {
	names(SI[[i]]) <- names(ast.grns[[i]])
}

## Calculate the % paralogs in each GRN
### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("../mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("../mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())]#[3:8]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))

names(go_lists) <- list.files()[grep("GO", list.files())]#[3:8]
names(go_lists)


paralog.numbers <- list()
paralog.numbers <- lapply(seq_along(conserved.grns), function(y) lapply(seq_along(conserved.grns[[y]]), function(x) calcParalog(conserved = conserved.grns[[y]], species.1 = zeb.grns[[y]], species.2 = ast.grns[[y]], ngenes.1 = length(go_lists[[13]][go_lists[[13]] %in% row.names(GetAssayData(hypo.zeb))]), ngenes.2 = length(go_lists[[13]][go_lists[[13]] %in% row.names(GetAssayData(hypo.ast))]), i = x)))

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

correlations <- lapply(seq_along(GO_lists), function(y) unlist(lapply(seq_along(GO_lists[[y]]), function(x) cor(GO_data[[y]][[x]]$weight_ast, GO_data[[y]][[x]]$weight_zeb))))

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

## Make SI plot for TFs vs NPs/NTS

SI.list <- readRDS("../SI_results_trinarized_a1.5_b2_f0.1.rds")
SI.sub.GO <- as.data.frame(SI.list$SI.sub.GO)
SI.plot <- ggplot(SI.sub.GO[SI.sub.GO$variable == c("NP_NTS", "TFs"),], aes(x = variable, y = value, group = variable, color = variable)) + geom_jitter(size = 1) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + xlab(element_blank()) + ylab("Similarity Index (SI)") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means()
# SI.plot <- ggplot(SI.sub.GO[grep("Neuronal_", SI.sub.GO$cell_type),], aes(x = variable, y = value, group = variable, color = variable)) + geom_jitter(size = 1) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + xlab(element_blank()) + ylab("Similarity Index (SI)") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means()
# SI.plot <- ggplot(SI.sub.GO, aes(x = variable, y = value, group = variable, color = variable)) + geom_jitter(size = 1) + geom_boxplot(outlier.color = NA, fill = "transparent", color = "black") + scale_color_viridis_d() + xlab(element_blank()) + ylab("Similarity Index (SI)") + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), axis.title = element_text(size = 10), legend.position = "none") + stat_compare_means()

SI.sub.GO.Neuronal_01 <- SI.list[[5]][grep("Neuronal_04", SI.list[[5]]$cell_type),]
SI.sub.GO.Neuronal_01$cell_type <- factor(SI.sub.GO.Neuronal_01$cell_type, levels = paste("Neuronal_04", c(0:11), sep = "_"))
SI.plot.Neuronal_01 <- ggplot(SI.sub.GO.Neuronal_01[13:36,], aes(x = cell_type, y = value, group = variable, color = variable)) + geom_point() + geom_line() + scale_color_viridis_d() + ylim(c(0,0.85)) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.title.x = element_blank(), legend.position = "none") + ylab("Similarity Index (SI)")

## Gather plots

# Make individual dot plots per gene

dot.plots <- GO_plots[["nps"]][["oxt"]] + GO_plots[["nts"]][["syt1a"]] + GO_plots[["synaptic"]][["slc32a1"]] + GO_plots[["ion"]][["kcnc1a"]] + GO_plots[["nps"]][["nmba"]] + GO_plots[["nts"]][["gna12a"]] + GO_plots[["synaptic"]][["syt1b"]] + GO_plots[["ion"]][["ryr1a"]] + plot_layout(nrow = 2, ncol = 4)

dot.plots <- dot.plots & theme(legend.position = "bottom")
dot.plots <- dot.plots + plot_layout(guides = "collect", nrow = 2)
dot.plots <- dot.plots & theme(axis.title = element_blank(), axis.text = element_text(size = 8), line = element_line(size = 1)) 

## Combine using patchwork, into something 8 inches wide by 3 inches tall
## NEED TO REMAKE THIS FOR CURRENT FIG LAYOUT

vip <- GO_plots[["nps"]][["vip"]] + theme(axis.text = element_text(size = 8), axis.title = element_blank(), legend.position = "none")



# Prep figures with proper sizes
fig_layout <- "
AB
CD
EF"

dev.new()
vip + scatter2 + plot_layout(width = unit(c(30,20), "mm"), height = unit(c(30), "mm"))
dev.new()
scatter1 + odds.scatter + plot_layout(width = unit(c(25,25), "mm"), height = unit(c(30), "mm"))
dev.new()
SI.plot + SI.plot.Neuronal_01 + plot_layout(width = unit(c(20,30), "mm"), height = unit(c(30), "mm"))

## Figure S6
dev.new()
dot.plots + plot_layout(heights = unit(c(35), c("mm")), widths = unit(c(35), c("mm")))

# 
# vip <- vip + plot_layout(width = unit(c(40), c("mm")), height = unit(c(40), c("mm")))
# scatter1 <- scatter1 + plot_layout(width = unit(c(28), c("mm")), height = unit(c(40), c("mm")))
# scatter2 <- scatter2 + plot_layout(width = unit(c(28), c("mm")), height = unit(c(40), c("mm")))
# 
# dev.new()
# vip + scatter1 + scatter2
# 
# 
# odds.scatter <- odds.scatter + plot_layout(width = unit(c(34), c("mm")), height = unit(c(40), c("mm")))
# SI.plot <- SI.plot + plot_layout(width = unit(c(28), c("mm")), height = unit(c(40), c("mm")))
# SI.plot.GABA_1 <- SI.plot.GABA_1 + plot_layout(width = unit(c(34), c("mm")), height = unit(c(40), c("mm")))
# 
# 
# dev.new()
# odds.scatter + SI.plot + SI.plot.GABA_1

