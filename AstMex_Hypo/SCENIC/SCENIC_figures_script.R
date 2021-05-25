library(reshape2)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(scales)
library(tictoc)
library(ggrepel)
library(viridis)
library(ggsignif)
library(ggpubr)
library(patchwork)

## I think I need to rewrite this and the functions script to account for more then two datasets, and datasets that aren't name ast and zeb (or surface and cave)
## Can use lapply on nested lists (by dataset, then by GO term, then by gene name)
## Adds a level of structure to this, but should make it simplier and easier
## Unclear how I will make plots, or which ones to make?


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")

# Load GENIE3 files from all the comps, for ast and zeb
datasets <- c("surface", "cave", "pachon", "tinaja", "molino")

gene.lists <- c("nps", "nts")

dataset.lists <- lapply(datasets, function(x) lapply(gene.lists, function(y) readRDS(paste(x, "/", y, "/int/1.4_GENIE3_linkList.Rds", sep = ""))))
names(dataset.lists) <- datasets

for (i in 1:length(dataset.lists)) {
  names(dataset.lists[[i]]) <- gene.lists
}

dataset.modules <- lapply(datasets, function(x) lapply(gene.lists, function(y) readRDS(paste(x, "/", y, "/int/1.6_tfModules_asDF.Rds", sep = ""))))
names(dataset.modules) <- datasets

for (i in 1:length(dataset.modules)) {
  names(dataset.modules[[i]]) <- gene.lists
}

# surface.lists <- lapply(gene.lists, function(x) readRDS(paste("surface/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
# cave.lists <- lapply(gene.lists, function(x) readRDS(paste("cave/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
# surface.modules <- lapply(gene.lists, function(x) readRDS(paste("surface/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
# cave.modules <- lapply(gene.lists, function(x) readRDS(paste("cave/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
# 
# pachon.lists <- lapply(gene.lists, function(x) readRDS(paste("pachon/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
# tinaja.lists <- lapply(gene.lists, function(x) readRDS(paste("tinaja/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
# molino.lists <- lapply(gene.lists, function(x) readRDS(paste("molino/", x, "/int/1.4_GENIE3_linkList.Rds", sep = "")))
# pachon.modules <- lapply(gene.lists, function(x) readRDS(paste("pachon/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
# tinaja.modules <- lapply(gene.lists, function(x) readRDS(paste("tinaja/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
# molino.modules <- lapply(gene.lists, function(x) readRDS(paste("molino/", x, "/int/1.6_tfModules_asDF.Rds", sep = "")))
# 
# names(surface.lists) <- gene.lists
# names(cave.lists) <- gene.lists
# names(surface.modules) <- gene.lists
# names(cave.modules) <- gene.lists
# 
# names(pachon.lists) <- gene.lists
# names(tinaja.lists) <- gene.lists
# names(molino.lists) <- gene.lists
# names(pachon.modules) <- gene.lists
# names(tinaja.modules) <- gene.lists
# names(molino.modules) <- gene.lists

## Make functions for extracting/plotting data
## load functions
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/SCENIC_functions.R")

## Load GO_lists and make plot data and plots

# extract the names of the genes belonging to each gene.list (results in a list of vectors length equal to gene.lists)
GO_lists <- sapply(gene.lists, function(x) Reduce(union, lapply(datasets, function(y) unique(dataset.lists[[y]][[x]]$Target))), USE.NAMES = T)

# Need to extract data for plotting
GO_data <- lapply(gene.lists, function(x) lapply(GO_lists[[x]], function(y) extractPlotData(linkLists = dataset.lists, modules = dataset.modules, go.term = x, gene.name = y)))

names(GO_data) <- names(GO_lists)
for (i in 1:length(GO_data)) {
  names(GO_data[[i]]) <- GO_lists[[i]]
}

# Make ggplots for each Target gene (terminal effector)

plotPlotData(data = GO_data$nps, gene = "hcrt", x = "surface", y = "cave", quantile = 0.95, log = F)

GO_plots_sc <- lapply(seq_along(GO_lists), function(x) lapply(names(GO_data[[x]]), function(y) plotPlotData(data = GO_data[[x]], gene = y, x = "surface", y = "cave", quantile = 0.95, log = F)))
GO_plots_sp <- lapply(seq_along(GO_lists), function(x) lapply(names(GO_data[[x]]), function(y) plotPlotData(data = GO_data[[x]], gene = y, x = "surface", y = "pachon", quantile = 0.95, log = F)))
GO_plots_st <- lapply(seq_along(GO_lists), function(x) lapply(names(GO_data[[x]]), function(y) plotPlotData(data = GO_data[[x]], gene = y, x = "surface", y = "tinaja", quantile = 0.95, log = F)))
GO_plots_sm <- lapply(seq_along(GO_lists), function(x) lapply(names(GO_data[[x]]), function(y) plotPlotData(data = GO_data[[x]], gene = y, x = "surface", y = "molino", quantile = 0.95, log = F)))

names(GO_plots_sc) <- names(GO_lists)
names(GO_plots_sp) <- names(GO_lists)
names(GO_plots_st) <- names(GO_lists)
names(GO_plots_sm) <- names(GO_lists)

for (i in 1:length(GO_plots_sc)) {
  names(GO_plots_sc[[i]]) <- GO_lists[[i]]
}
for (i in 1:length(GO_plots_sp)) {
	names(GO_plots_sp[[i]]) <- GO_lists[[i]]
}
for (i in 1:length(GO_plots_st)) {
  names(GO_plots_st[[i]]) <- GO_lists[[i]]
}
for (i in 1:length(GO_plots_sm)) {
  names(GO_plots_sm[[i]]) <- GO_lists[[i]]
}

sc <- GO_plots_sc[["nps"]][["hcrt"]]
sp <- GO_plots_sp[["nps"]][["hcrt"]]
st <- GO_plots_st[["nps"]][["hcrt"]]
sm <- GO_plots_sm[["nps"]][["hcrt"]]

plots <- sc + sp + st + sm + plot_layout(ncol = 4)
plots <- plots & xlim(c(-0.04, 0.04))
plots <- plots & ylim(c(-0.04, 0.04))
plots

## Calculate the Drift Index for each NP for plotting
## Use weight cutoffs for calling GRNs, 0.001 and 0.005 or whether the target is in the top50 for each TF (top50), or the top50 TFs per target (top50perTarget)
## So tfModules_asDF.zeb and tfModules_asDF.ast have done this for me already

surface.grns <- lapply(surface.modules, function(x) x[x$method == "top50",])

for (i in 1:length(surface.grns)) {
	names <- unique(surface.grns[[i]]$Target)
	surface.grns[[i]] <- lapply(names, function(x) surface.grns[[i]][surface.grns[[i]]$Target == x, "TF"])
	names(surface.grns[[i]]) <- names
}

mean(unlist(lapply(cave.grns, function(x) lapply(x, function(y) length(y)))))

cave.grns <- lapply(cave.modules, function(x) x[x$method == "top50",])

for (i in 1:length(cave.grns)) {
	names <- unique(cave.grns[[i]]$Target)
	cave.grns[[i]] <- lapply(names, function(x) cave.grns[[i]][cave.grns[[i]]$Target == x, "TF"])
	names(cave.grns[[i]]) <- names
	cave.grns <- cave.grns[intersect(names(surface.grns), names(cave.grns))]
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


## Calculate the Drift Index

DI <- lapply(seq_along(conserved.grns), function(x) calcDriftIndex(conserved = conserved.grns[[x]], species.1 = surface.grns[[x]], species.2 = cave.grns[[x]]))

names(DI) <- names(cave.grns)
for (i in 1:length(DI)) {
	names(DI[[i]]) <- names(cave.grns[[i]])
}

## Calculate the % paralogs in each GRN
### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Otophysi_Paralogs/mart_export_Danio.txt", sep = "\t", head = TRUE)
mart[[2]] <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Otophysi_Paralogs/mart_export_Astyanax_102.txt", sep = "\t", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())][3:8]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))

names(go_lists) <- list.files()[grep("GO", list.files())][3:8]
names(go_lists)


paralog.numbers <- list()
paralog.numbers <- lapply(seq_along(conserved.grns), function(y) lapply(seq_along(conserved.grns[[y]]), function(x) calcParalog(conserved = conserved.grns[[y]], species.1 = surface.grns[[y]], species.2 = cave.grns[[y]], ngenes.1 = length(go_lists[[1]][go_lists[[1]] %in% row.names(GetAssayData(hypo.ast))]), ngenes.2 = length(go_lists[[1]][go_lists[[1]] %in% row.names(GetAssayData(hypo.ast))]), i = x)))

grn.paralog <- lapply(paralog.numbers, function(x) data.frame(do.call(rbind, x)))
grn.paralog <- lapply(grn.paralog, function(x) rbind(x, colSums(x)))

for (i in 1:length(grn.paralog)) {
	row.names(grn.paralog[[i]]) <- c(names(conserved.grns[[i]]), "sums")
	grn.paralog[[i]]$percent.para <- unlist(apply(grn.paralog[[i]], 1, function(x) sum(x[1], x[5])/sum(x[2], x[6]))*100)
	grn.paralog[[i]]$percent.conserved <- unlist(apply(grn.paralog[[i]], 1, function(x) x[9]/x[10]*100))
	# grn.paralog[[i]]$fisher.pval <- as.numeric(fisher.results[,1])
	grn.paralog[[i]]$fisher.odds <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[1:4] + x[5:8], ncol = 2, nrow = 2))))[,3])
	grn.paralog[[i]]$fisher.odds.pval <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[1:4] + x[5:8], ncol = 2, nrow = 2))))[,1])
	grn.paralog[[i]]$fisher.odds.2 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2))))[,3])
	grn.paralog[[i]]$fisher.odds.pval.2 <- as.numeric(do.call(rbind, apply(grn.paralog[[i]], 1, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2))))[,1])
	grn.paralog[[i]]$gene <- row.names(grn.paralog[[i]])
}


## Calculate the correlation between the GENIE3 results for each species, for each GOI

correlations <- lapply(seq_along(GO_lists), function(y) unlist(lapply(seq_along(GO_lists[[y]]), function(x) cor(GO_data[[y]][[x]]$weight_ast, GO_data[[y]][[x]]$weight_zeb))))

names(correlations) <- names(GO_lists)
for (i in 1:length(correlations)) {
	names(correlations[[i]]) <- GO_lists[[i]]
}

## Put Correlation and DI into the same df
## gene / GO / DI / corr

DI <- lapply(DI, function(x) unlist(x))
correlations <- lapply(correlations, function(x) unlist(x))
correlations <- lapply(seq_along(correlations), function(x) correlations[[x]][names(DI[[x]])])
names(correlations) <- names(DI)


df <- melt(do.call(rbind, lapply(seq_along(DI), function(x) data.frame(gene = names(DI[[x]]), GO = names(DI)[x], DI = DI[[x]], corr = correlations[[x]], percent.para = grn.paralog[[x]]$percent.para[1:length(DI[[x]])]))))

df2 <- do.call(rbind, lapply(seq_along(DI), function(x) data.frame(gene = names(DI[[x]]), GO = names(DI)[x], fisher.odds = grn.paralog[[x]]$fisher.odds[1:length(DI[[x]])], fisher.odds.pval = grn.paralog[[x]]$fisher.odds.pval[1:length(DI[[x]])])))

# Re-order factor levels
# Some genes are in multiple categories, and therefore can't be factored correctly...

df <- df[!duplicated(df[c(1,3)]),]

df$gene <- factor(df$gene, levels = df[order(df$value[df$variable == "corr"]), "gene"], ordered = T)
df$variable <- factor(df$variable, levels = c("corr", "DI", "percent.para"))
## Plot scatter plots and geom_tile heatmaps!

anno_df <- compare_means(value ~ GO, group.by = 'variable', data = df)
anno_df$y_pos <- c( 0.75, 0.9, 0.95, 1.0, 1.05, 1.1, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 15, 16, 17, 18, 19, 20)

scatter <- ggplot(df, aes(x = GO, y = value, color = variable)) + geom_point(position = position_jitterdodge(), size = .5) + geom_boxplot(fill = NA, size = .5, outlier.shape = NA, color = "black") + scale_color_manual(values = c("grey55", "red", "black")) + theme(axis.title = element_blank(), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8)) + scale_x_discrete(labels = c("Neuropeptides", "Neurotransmitters", "Synaptic Genes", "Ion Channel")) + facet_wrap(~variable, scales = "free") + guides(color = F) 
+ geom_signif(data = anno_df, aes(xmin = group1, xmax = group2, annotations = p.signif, y_position = y_pos), manual = T, textsize = 3)

scatter2 <- ggplot(df2, aes(x = GO, y = fisher.odds, color = ifelse(fisher.odds.pval < 0.05, "yes", "no"))) + geom_jitter(size = 0.5) + geom_boxplot(fill = NA, size = 0.5, outlier.shape = NA, color = "black") + scale_color_manual(values = c("black", "red")) + theme(axis.title = element_text(size = 10), legend.title = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8)) + scale_x_discrete(labels = c("Neuropeptides", "Neurotransmitters", "Synaptic Genes", "Ion Channel")) + ylab("Odds ratio") + xlab("") + guides(color = F)


# Make heatmaps

corr.plot <- ggplot(df[df$variable == "corr",], aes(x = gene, y = variable, fill = value)) + geom_tile() + facet_grid(~GO, space = "free", scales = "free") + theme(axis.text.x = element_blank(), axis.title.x = element_blank()) + scale_fill_viridis(option = "A")

di.plot <- ggplot(df[df$variable == "DI",], aes(x = gene, y = variable, fill = value)) + geom_tile() + facet_grid(~GO, space = "free", scales = "free") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) + scale_fill_viridis(option = "A")

heatmaps <- corr.plot / di.plot + plot_layout(guides = "collect")
heatmaps <- heatmaps & theme(axis.title = element_blank())

GO_plots[["nps"]][["hcrt"]]

# Make individual dot plots per gene

dot.plots <- GO_plots[["nps"]][["oxt"]] + GO_plots[["nts"]][["isl1"]] + GO_plots[["synaptic"]][["slc32a1"]] + GO_plots[["ion"]][["grin1b"]] + GO_plots[["nps"]][["nmba"]] + GO_plots[["nts"]][["pfas"]] + GO_plots[["synaptic"]][["snap25b"]] + GO_plots[["ion"]][["anxa6"]]

dot.plots <- dot.plots & theme(legend.position = "bottom")
dot.plots <- dot.plots + plot_layout(guides = "collect", nrow = 2)
dot.plots <- dot.plots & theme(axis.title = element_blank(), axis.text = element_text(size = 5)) 


## Combine using patchwork, into something 8 inches wide by 3 inches tall

patch <- 
((heatmaps[[1]]/heatmaps[[2]]) | scatter + plot_layout(widths = c(1,.5))) / dot.plots 

patch <- (scatter + dot.plots + plot_layout(widths = c(1,3))) / (heatmaps) + plot_layout(heights = unit(c(3,1.5), c('in', 'in')), widths = unit(c(7), c('in')))

patch[[2]] <- patch[[2]] + theme(axis.text.x = element_text(size = 4), legend.text = element_text(size = 8))

patch


