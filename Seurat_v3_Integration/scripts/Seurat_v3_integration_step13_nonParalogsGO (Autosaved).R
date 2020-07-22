library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(data.table)
library(fdrtool)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")


## Functions

splitNames <- function(x) {
	if (length(x) == 2) {
		names(x) <- c("danio", "astyanax")
	}
	return(x)
}

listList <- function(x) {
	list <- list()
	for (i in 1:length(x)){
		list[[i]] <- lapply(x[[i]], function(y) length(y) > 0)
	}
	return(list)
}

hCluster <- function(x = go_analysis, measure.var = "Benj.value", category = "KEGG_PATHWAY") {
	casted <- acast(x[x$Category == category,], formula = L1 + L2 ~ Term, value.var = measure.var)
	casted[is.na(casted)] <- 0
	hr <- hclust(as.dist(1-cor(t(casted), method="pearson")), method="complete")
	hc <- hclust(as.dist(1-cor(casted, method="spearman")), method="complete")
	clustered <- x[x$Category == category,]
	clustered$Term <- factor(clustered$Term, levels = hc[[4]][hc[[3]]])
	clustered$L3 <- factor(clustered$L3, levels = hr[[4]][hr[[3]]])
	return(clustered)
}

## Load gene lists

gene.lists <- readRDS(file = "drift_gene_lists.rds")

names(gene.lists) <- c("conserved.markers", "zebrafish.markers", "astyanax.markers", "conserved.markers.sub", "zebrafish.markers.sub", "astyanax.markers.sub", "conserved.markers.ast", "surface.markers", "cave.markers", "conserved.markers.ast.sub", "surface.markers.sub", "cave.markers.sub")

## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4)] <- lapply(gene.lists[c(1,4)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$zebrafish_avg_logFC > 0 & x[[y]]$astyanax_avg_logFC > 0,]))
gene.lists.pos[c(2,3,5,6,8,9,11,12)] <- lapply(gene.lists[c(2,3,5,6,8,9,11,12)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))
gene.lists.pos[c(7,10)] <- lapply(gene.lists[c(7,10)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}



### Make new gene list that is all the subcluster genes, but named by subtype

str(gene.lists.pos, max.level = 1)

# for loop, or apply over the names of subtypes, and do this for all three gene.lists.pos[[4:6]]

test <- lapply(list(4,5,6), function(y) lapply(names(gene.lists.pos[[1]]), function(x) Reduce(union, lapply(gene.lists.pos[[y]][grep(x, names(gene.lists.pos[[4]]))], function(z) row.names(z)))))
names(test) <- c("conserved.markers.sub.combined", "zebrafish.markers.sub.combined", "astyanax.markers.sub.combined")


for(i in 1:length(test)){
	names(test[[i]]) <- names(gene.lists[[i]])
}

gene.lists.pos <- c(gene.lists.pos, test)


## Take setdiff for species specific markers

annoCharts <- list()
for (i in 1:length(gene.lists.pos[[2]])) {
	markers <- unique(setdiff(gene.lists.pos[[14]][[i]], gene.lists.pos[[13]][[i]]))
	# markers <- markers[markers %in% row.names(hypo.zeb@data)]
	# markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(gene.lists.pos[[2]])[[i]], listType = "Gene")
	annoCharts[[i]] <- list()
	annoCharts[[i]][[1]] <- getFunctionalAnnotationChart(david)
	markers <- unique(setdiff(gene.lists.pos[[15]][[i]], gene.lists.pos[[13]][[i]]))
	# markers <- markers[markers %in% row.names(hypo.zeb@data)]
	# markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(gene.lists.pos[[3]])[i], listType = "Gene")
	annoCharts[[i]][[2]] <- getFunctionalAnnotationChart(david)
}

# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts) <- names(gene.lists.pos[[2]])
annoCharts <- lapply(annoCharts, function(x) splitNames(x))

list.index <- listList(annoCharts)
annoCharts.2 <- lapply(seq_along(annoCharts), function(x) annoCharts[[x]][unlist(list.index[[x]])])
names(annoCharts.2) <- names(gene.lists.pos[[2]])
annoCharts.2 <- lapply(annoCharts.2, function(x) splitNames(x))

# Reshape for plotting
measure.vars <- matrix(c("Bonferroni", "Benjamini", "FDR", "Count", "Fold.Enrichment", "X.", "Bonferroni", "Benjamini", "FDR", "Count", "FE", "X.", "Bonf.value", "Benj.value", "FDR.value", "Counts.value", "FoldE.value", "X.value"), nrow=6, ncol = 3)
melted <- apply(measure.vars, 1, function(x) melt(annoCharts.2, id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))

go_analysis <- melted[[1]]
go_analysis$Benj.value <- melted[[2]]$Benj.value
go_analysis$FDR.value <- melted[[3]]$FDR.value
go_analysis$Count <- melted[[4]]$Counts.value
go_analysis$FoldE <- log(melted[[5]]$FoldE.value)
go_analysis$Perc <- melted[[6]]$X.value
go_analysis$FoldE[go_analysis$L2 == "astyanax"] <- go_analysis$FoldE[go_analysis$L2 == "astyanax"]*-1
go_analysis$L3 <- paste(go_analysis$L1, go_analysis$L2, sep = "_")

# Cluster with function, and plot
# GOTERM_BP_DIRECT
# UP_KEYWORDS
# GOTERM_MF_DIRECT
# GOTERM_CC_DIRECT
# KEGG_PATHWAY

png("Figures/GO_analysis_nonparalogs_subtypes.png", width = 12, height = 6, units = "in", res = 250)
ggplot(data = hCluster(x = go_analysis[go_analysis$Benj.value < .05,], measure.var = "Benj.value", category = "UP_KEYWORDS"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point(alpha = 0.75) + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(size = 6), axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()

## Do same for subclusters

annoCharts.sub <- list()
for (i in 1:length(gene.lists.pos[[4]])) {
	markers <- unique(setdiff(row.names(gene.lists.pos[[5]][[i]]), row.names(gene.lists.pos[[4]][[i]])))
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(gene.lists.pos[[5]])[[i]], listType = "Gene")
	annoCharts.sub[[i]] <- list()
	annoCharts.sub[[i]][[1]] <- getFunctionalAnnotationChart(david)
	markers <- unique(setdiff(row.names(gene.lists.pos[[6]][[i]]), row.names(gene.lists.pos[[4]][[i]])))
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(gene.lists.pos[[6]])[i], listType = "Gene")
	annoCharts.sub[[i]][[2]] <- getFunctionalAnnotationChart(david)
}

# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts.sub) <- names(gene.lists.pos[[4]])
annoCharts.sub <- lapply(annoCharts.sub, function(x) splitNames(x))

list.index <- listList(annoCharts.sub)
annoCharts.sub.2 <- lapply(seq_along(annoCharts.sub), function(x) annoCharts.sub[[x]][unlist(list.index[[x]])])
names(annoCharts.sub.2) <- names(gene.lists.pos[[4]])
annoCharts.sub.2 <- lapply(annoCharts.sub.2, function(x) splitNames(x))

# Reshape for plotting
measure.vars <- matrix(c("Bonferroni", "Benjamini", "FDR", "Count", "Fold.Enrichment", "X.", "Bonferroni", "Benjamini", "FDR", "Count", "FE", "X.", "Bonf.value", "Benj.value", "FDR.value", "Counts.value", "FoldE.value", "X.value"), nrow=6, ncol = 3)
melted <- apply(measure.vars, 1, function(x) melt(annoCharts.sub.2, id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))

go_analysis.sub <- melted[[1]]
go_analysis.sub$Benj.value <- melted[[2]]$Benj.value
go_analysis.sub$FDR.value <- melted[[3]]$FDR.value
go_analysis.sub$Count <- melted[[4]]$Counts.value
go_analysis.sub$FoldE <- log(melted[[5]]$FoldE.value)
go_analysis.sub$Perc <- melted[[6]]$X.value
go_analysis.sub$FoldE[go_analysis.sub$L2 == "astyanax"] <- go_analysis.sub$FoldE[go_analysis.sub$L2 == "astyanax"]*-1
go_analysis.sub$L3 <- paste(go_analysis.sub$L1, go_analysis.sub$L2, sep = "_")
go_analysis.sub$L4 <- hypo.integrated@meta.data$integrated_Subtype[match(go_analysis.sub$L1, hypo.integrated@meta.data$integrated_SubclusterType)]
go_analysis.sub$L5 <- paste(go_analysis.sub$L4, go_analysis.sub$L2, sep = "_")

# Cluster with function, and plot
# GOTERM_BP_DIRECT
# UP_KEYWORDS
# GOTERM_MF_DIRECT
# GOTERM_CC_DIRECT 
# KEGG_PATHWAY

ggplot(data = hCluster(x = go_analysis.sub[go_analysis.sub$FDR.value < 1,], measure.var = "Benj.value", category = "UP_KEYWORDS"), aes(L4, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point(alpha = 0.5) + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text = element_text(size = 8), axis.title.y = element_blank()) + xlab("UP_KEYWORDS Term") + coord_flip()


### Save annoCharts for later, even though they don't show much!

anno.list <- c(annoCharts.2, annoCharts.sub.2)
saveRDS(anno.list, file = "GO_analysis_all.rds")
