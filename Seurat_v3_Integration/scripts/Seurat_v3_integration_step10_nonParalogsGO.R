library(Matrix)
library(dplyr)
library(ggplot2)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(data.table)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

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

### Make new gene list that is all the subcluster genes, but named by cluster
gene.lists.pos <- readRDS("drift_gene_lists_pos_trinarized_a1.5_b2_f0.1.rds")

## Subset this list for non-paralogs

# Need to do this for (1) 2,3, and (7) 8,9

gene.lists.non.para <- gene.lists.pos

clusters <- names(gene.lists.pos[[1]])
subclusters <- names(gene.lists.pos[[7]])

gene.lists.non.para[[2]] <- lapply(seq_along(clusters), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[2]][[clusters[x]]]), row.names(gene.lists.pos[[1]][[clusters[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[3]][[clusters[x]]]), row.names(gene.lists.pos[[1]][[clusters[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[clusters[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[clusters[x]]]), mart[[2]]$Gene.name)])
  paralog.2 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart[[2]]$Gene.name)])
  paralog.2 <- union(paralog.con, paralog.2)
  paralog.2 <- paralog.2[paralog.2 %in% mart[[1]]$Gene.name]
  return(gene.lists.pos[[2]][[clusters[x]]][genes.1[!(genes.1 %in% paralog.2)],])
})
gene.lists.non.para[[3]] <- lapply(seq_along(clusters), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[2]][[clusters[x]]]), row.names(gene.lists.pos[[1]][[clusters[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[3]][[clusters[x]]]), row.names(gene.lists.pos[[1]][[clusters[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[clusters[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[clusters[x]]]), mart[[2]]$Gene.name)])
  paralog.1 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart[[2]]$Gene.name)])
  paralog.1 <- union(paralog.con, paralog.1)
  paralog.1 <- paralog.1[paralog.1 %in% mart[[2]]$Gene.name]
  return(gene.lists.pos[[3]][[clusters[x]]][genes.2[!(genes.2 %in% paralog.1)],])
})

names(gene.lists.non.para[[2]]) <- clusters
names(gene.lists.non.para[[2]]) <- clusters

gene.lists.non.para[[8]] <- lapply(seq_along(subclusters), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[8]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[9]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[2]]$Gene.name)])
  paralog.2 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart[[2]]$Gene.name)])
  paralog.2 <- union(paralog.con, paralog.2)
  paralog.2 <- paralog.2[paralog.2 %in% mart[[1]]$Gene.name]
  return(gene.lists.pos[[8]][[subclusters[x]]][genes.1[!(genes.1 %in% paralog.2)],])
})
gene.lists.non.para[[9]] <- lapply(seq_along(subclusters), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[8]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[9]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[2]]$Gene.name)])
  paralog.1 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart[[2]]$Gene.name)])
  paralog.1 <- union(paralog.con, paralog.1)
  paralog.1 <- paralog.1[paralog.1 %in% mart[[2]]$Gene.name]
  return(gene.lists.pos[[9]][[subclusters[x]]][genes.2[!(genes.2 %in% paralog.1)],])
})

names(gene.lists.non.para[[8]]) <- subclusters
names(gene.lists.non.para[[9]]) <- subclusters

# for loop, or apply over the names of clusters, and do this for all three gene.lists.pos[[4:6]]
test <- lapply(list(7,8,9), function(y) lapply(names(gene.lists.non.para[[1]]), function(x) Reduce(union, lapply(gene.lists.non.para[[y]][grep(x, names(gene.lists.non.para[[y]]))], function(z) row.names(z)))))
names(test) <- c("conserved.markers.sub.combined", "zebrafish.markers.sub.combined", "astyanax.markers.sub.combined")


for(i in 1:length(test)){
	names(test[[i]]) <- names(gene.lists.pos[[1]])
}

## Append gene lists
gene.lists.pos <- c(gene.lists.pos[1:14], test)


## Take setdiff for species specific markers

annoCharts <- list()
for (i in 1:length(gene.lists.pos[["conserved.markers.sub.combined"]])) {
	markers <- unique(setdiff(gene.lists.pos[["zebrafish.markers.sub.combined"]][[i]], gene.lists.pos[["conserved.markers.sub.combined"]][[i]]))
	# markers <- markers[markers %in% row.names(hypo.zeb@data)]
	# markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "BBID"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(gene.lists.pos[[1]])[[i]], listType = "Gene")
	annoCharts[[i]] <- list()
	annoCharts[[i]][[1]] <- getFunctionalAnnotationChart(david)
	markers <- unique(setdiff(gene.lists.pos[["astyanax.markers.sub.combined"]][[i]], gene.lists.pos[["conserved.markers.sub.combined"]][[i]]))
	# markers <- markers[markers %in% row.names(hypo.zeb@data)]
	# markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "UP_SEQ_FEATURE", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY", "REACTOME_PATHWAY", "BBID"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(gene.lists.pos[[1]])[i], listType = "Gene")
	annoCharts[[i]][[2]] <- getFunctionalAnnotationChart(david)
}

# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts) <- names(gene.lists.pos[["conserved.markers.sub.combined"]])
annoCharts <- lapply(annoCharts, function(x) splitNames(x))

list.index <- listList(annoCharts)
annoCharts.2 <- lapply(seq_along(annoCharts), function(x) annoCharts[[x]][unlist(list.index[[x]])])
names(annoCharts.2) <- names(gene.lists.pos[["conserved.markers.sub.combined"]])
annoCharts.2 <- lapply(annoCharts.2, function(x) splitNames(x))

saveRDS(annoCharts.2, file = "GO_analysis_all_trinarized.rds")

# # If used subcluster markers, but grouped by cluster
annoCharts.2 <- readRDS("GO_analysis_all_trinarized.rds")

# Reshape for plotting
measure.vars <- matrix(c("Bonferroni", "Benjamini", "FDR", "Count", "Fold.Enrichment", "X.", "Bonferroni", "Benjamini", "FDR", "Count", "FE", "X.", "Bonf.value", "Benj.value", "FDR.value", "Counts.value", "FoldE.value", "X.value"), nrow=6, ncol = 3)
melted <- apply(measure.vars, 1, function(x) reshape2::melt(annoCharts.2, id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))

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
# UP_SEQ_FEATURE
# REACTOME_PATHWAY
# BBID

go_analysis$L1 <- factor(go_analysis$L1, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

go.plot <- ggplot(data = hCluster(x = go_analysis[go_analysis$Benj.value < .05,], measure.var = "Benj.value", category = "KEGG_PATHWAY"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point(alpha = 0.75) + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_blank(), axis.title.y = element_blank()) + coord_flip()

pdf("Figures/Hypo_integrated_GO-enrichment_plot.pdf", height = 11, width = 9) 
go.plot + plot_layout(width = unit(120, "mm"), height = unit(90, "mm"))
dev.off()


