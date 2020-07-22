library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)

# Load subsets
load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj")
load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_65k.Robj")

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")
markers.Subtype.species <- readRDS("Hypo_integrated_markers.species.SubType.list.rds")
markers.SubclusterType.species <- readRDS("Hypo_integrated_markers.species.SubclusterType.list.rds")

# # Create the david object, associated with registered email address
# david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")

# # Examine object, and any gene lists associated with it
# david

# getIdTypes(david)
# getAllAnnotationCategoryNames(david)

splitMarkers <- function(x) {
	x.list <- list()
	x.list[[1]] <- x[x$avg_logFC > 1,]
	x.list[[2]] <- x[x$avg_logFC < -1,]
	return(x.list)
}

splitNames <- function(x) {
	names(x) <- c("danio", "astyanax")
	return(x)
}

markers.Subtype.species <- lapply(markers.Subtype.species, function(x) splitMarkers(x))

annoCharts <- list()
for (i in 1:length(markers.Subtype.species)) {
	markers <- row.names(markers.Subtype.species[[i]][[1]])
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.Subtype.species)[[i]], listType = "Gene")
	annoCharts[[i]] <- list()
	annoCharts[[i]][[1]] <- getFunctionalAnnotationChart(david)
	markers <- row.names(markers.Subtype.species[[i]][[2]])
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.Subtype.species)[i], listType = "Gene")
	annoCharts[[i]][[2]] <- getFunctionalAnnotationChart(david)
}

# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts) <- names(markers.Subtype.species)
annoCharts <- lapply(annoCharts, function(x) splitNames(x))

saveRDS(annoCharts, file = "GO_analysis_Subtypes.rds")


# Same for subcluster markers

markers.SubclusterType.species <- lapply(markers.SubclusterType.species[setdiff(1:length(markers.SubclusterType.species), grep("species specific cell type", markers.SubclusterType.species))], function(x) splitMarkers(as.data.frame(x)))

annoCharts.sc <- list()
for (i in 47:length(markers.SubclusterType.species)) {
	markers <- row.names(markers.SubclusterType.species[[i]][[1]])
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	if(length(markers) > 0) {
		markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
		david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
		setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
		result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.SubclusterType.species)[i], listType = "Gene")
		annoCharts.sc[[i]] <- list()
		annoCharts.sc[[i]][[1]] <- getFunctionalAnnotationChart(david)
	} else {
	annoCharts.sc[[i]] <- list()
	annoCharts.sc[[i]][[1]] <- "no signficantly differentially expressed genes"
	}
	markers <- row.names(markers.SubclusterType.species[[i]][[2]])
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	if(length(markers) > 0) {
		markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
		david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
		setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
		result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.SubclusterType.species)[i], listType = "Gene")
		annoCharts.sc[[i]][[2]] <- getFunctionalAnnotationChart(david)
	} else {
	annoCharts.sc[[i]] <- list()
	annoCharts.sc[[i]][[2]] <- "no signficantly differentially expressed genes"
	}
}

# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts.sc) <- names(markers.SubclusterType.species)

annoCharts.sc <- lapply(annoCharts.sc, function(x) splitNames(x))

saveRDS(annoCharts.sc, file = "GO_analysis_SubclusterTypes.rds")


# Reshape for plotting
measure.vars <- matrix(c("Bonferroni", "Benjamini", "FDR", "Count", "Fold.Enrichment", "X.", "Bonferroni", "Benjamini", "FDR", "Count", "FE", "X.", "Bonf.value", "Benj.value", "FDR.value", "Counts.value", "FoldE.value", "X.value"), nrow=6, ncol = 3)
melted <- apply(measure.vars, 1, function(x) melt(annoCharts, id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))

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
# GOTERM_CC_DIRECT KEGG_PATHWAY

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

ggplot(data = hCluster(x = go_analysis, measure.var = "Benj.value", category = "GOTERM_CC_DIRECT"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point() + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))

ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes_UP_KEYWORDS.png", units = "in", height = 10, width = 10, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes_KEGG_PATHWAY.png", units = "in", height = 10, width = 10, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes_GOTERM_BP_DIRECT.png", units = "in", height = 20, width = 18, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes_GOTERM_MF_DIRECT.png", units = "in", height = 15, width = 18, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes_GOTERM_CC_DIRECT.png", units = "in", height = 10, width = 12, dpi = 250)



index <- vector()
for (i in 1:length(markers.SubclusterType.species)) {
	markers <- row.names(markers.SubclusterType.species[[i]][[1]])
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	if(length(markers) > 0) {
		index[i] <- TRUE
	} else {
		index[i] <- FALSE
	}
	markers <- row.names(markers.SubclusterType.species[[i]][[2]])
	markers <- markers[markers %in% row.names(hypo.zeb@data)]
	markers <- markers[markers %in% row.names(hypo.ast@data)]
	if(length(markers) > 0) {
		index[i] <- TRUE
	} else {
		index[i] <- FALSE
	}
}


# Subclusters
# Reshape for plotting
melted <- apply(measure.vars, 1, function(x) melt(annoCharts.sc[index], id.vars = c("Category", "Term"), measure.vars = x[1], variable.name = x[2], value.name = x[3]))

go_analysis.sc <- melted[[1]]
go_analysis.sc$Benj.value <- melted[[2]]$Benj.value
go_analysis.sc$FDR.value <- melted[[3]]$FDR.value
go_analysis.sc$Count <- melted[[4]]$Counts.value
go_analysis.sc$FoldE <- log(melted[[5]]$FoldE.value)
go_analysis.sc$Perc <- melted[[6]]$X.value
go_analysis.sc$FoldE[go_analysis.sc$L2 == "astyanax"] <- go_analysis.sc$FoldE[go_analysis.sc$L2 == "astyanax"]*-1
go_analysis.sc$L3 <- paste(go_analysis.sc$L1, go_analysis.sc$L2, sep = "_")


ggplot(data = hCluster(x = go_analysis.sc, measure.var = "Benj.value", category = "GOTERM_CC_DIRECT"), aes(L1, Term, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point() + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))

ggsave("Figures/Zeb_Ast_DE_GO_analysis_SubclusterTypes_UP_KEYWORDS.png", units = "in", height = 20, width = 35, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_SubclusterTypes_KEGG_PATHWAY.png", units = "in", height = 10, width = 35, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_SubclusterTypes_GOTERM_BP_DIRECT.png", units = "in", height = 40, width = 35, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_SubclusterTypes_GOTERM_MF_DIRECT.png", units = "in", height = 40, width = 35, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_SubclusterTypes_GOTERM_CC_DIRECT.png", units = "in", height = 20, width = 35, dpi = 250)



