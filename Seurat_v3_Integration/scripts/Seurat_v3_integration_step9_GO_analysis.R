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
setwd("~/Documents/Schier_Lab/R_Projects/Seurat_v3_Integration/")
markers.Subtypes.filtered <- readRDS("Zeb_Ast_markers.Subtypes.filtered.rds")
markers.Subtypes.filtered <- readRDS("Zeb_Ast_markers.SubclusterTypes.filtered.rds")

# Create the david object, associated with registered email address
david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")

# Examine object, and any gene lists associated with it
david

getIdTypes(david)
getAllAnnotationCategoryNames(david)

splitNames <- function(x) {
	names(x) <- c("danio", "astyanax")
	return(x)
}

markers.Subtypes.filtered.species <- lapply(markers.Subtypes.filtered, function(x) splitMarkers(x))

annoCharts <- list()
for (i in 1:length(markers.Subtypes.filtered.species)) {
	markers <- row.names(markers.Subtypes.filtered.species[[i]][[1]])
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.Subtypes.filtered.species)[i], listType = "Gene")
	annoCharts[[i]] <- list()
	annoCharts[[i]][[1]] <- getFunctionalAnnotationChart(david)
	markers <- row.names(markers.Subtypes.filtered.species[[i]][[2]])
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.Subtypes.filtered.species)[i], listType = "Gene")
	annoCharts[[i]][[2]] <- getFunctionalAnnotationChart(david)
}

# Above makes a list of list, melt considers each list, so is Null? Maybe lapply melt, then melt the result?

names(annoCharts) <- names(markers.Subtypes.filtered)

annoCharts <- lapply(annoCharts, function(x) splitNames(x))

Benj <- melt(annoCharts, id.vars = c("Category", "Term"), measure.vars = c("Benjamini"), variable.name = "Benjamini", value.name = "Benj.value")
Count <- melt(annoCharts, id.vars = c("Category", "Term"), measure.vars = c("Count"), variable.name = "Count", value.name = "Counts.value")
FoldE <- melt(annoCharts, id.vars = c("Category", "Term"), measure.vars = c("Fold.Enrichment"), variable.name = "FE", value.name = "FoldE.value")
Perc <- melt(annoCharts, id.vars = c("Category", "Term"), measure.vars = c("X."), variable.name = "X.", value.name = "X.value")

test2 <- Benj
test2$Count <- Count$Counts.value
test2$FoldE <- log(FoldE$FoldE.value)
test2$Perc <- Perc$X.value


test2$FoldE[test2$L2 == "astyanax"] <- test2$FoldE[test2$L2 == "astyanax"]*-1
test2$L3 <- paste(test2$L1, test2$L2, sep = "_")


# cast it back to do clustering?
# Cluster based on Benj value, or whatever else
# Need to modify this so that I'm clustering only those terms from a given category, and not all of them?
# Fix to add Category with Term??

test <- acast(Benj, formula = L1 + L2 ~ Term, value.var="Benj.value")
test[is.na(test)] <- 0

hr <- hclust(as.dist(1-cor(t(test), method="pearson")), method="complete")

hc <- hclust(as.dist(1-cor(test, method="spearman")), method="complete")

# hr and hc contain the correct order in the 3rd list, of the names in the 4th list

L1.order <- hr[[4]][hr[[3]]]
Term.order <- hc[[4]][hc[[3]]]


# Plot as before, but with new row and column orders
test2$Term <- factor(test2$Term, levels = Term.order)
test2$L3 <- factor(test2$L3, levels = L1.order)

ggplot(data = test2[test2$Category == "UP_KEYWORDS",], aes(Term, L1, size = Benj.value, fill = FoldE, shape = L2, color = L2)) + geom_point() + scale_shape_manual(values = c(21, 24)) + scale_size_continuous(range = c(6, 2)) + scale_color_manual(values = c("black", "black")) + scale_fill_gradient2(low = "skyblue2", high = "khaki2") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 4, hjust = 1))



ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes.png", units = "in", height = 30, width = 10, dpi = 250)
ggsave("Figures/Zeb_Ast_DE_GO_analysis_Subtypes_UP_KEYWORDS.png", units = "in", height = 15, width = 7.5, dpi = 250)

ggsave("Figures/Zeb_Ast_DE_GO_analysis_SubclusterTypes.png", units = "in", height = 25, width = 25, dpi = 250)








annoCharts <- list()
for (i in setdiff(1:length(markers.SubclusterTypes.filtered), grep("Too few cells for comparison of", markers.SubclusterTypes.filtered))) {
	markers <- row.names(markers.SubclusterTypes.filtered[[i]])
	markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
	david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
	setAnnotationCategories(david, c("KEGG_PATHWAY"))
	result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = names(markers.SubclusterTypes.filtered)[i], listType = "Gene")
	annoCharts[[i]] <- getFunctionalAnnotationChart(david)
}

# The list will have NULL entries for missing SubclusterTypes, remove them
annoCharts <- annoCharts[c(setdiff(1:length(markers.SubclusterTypes.filtered), grep("Too few cells for comparison of", markers.SubclusterTypes.filtered)))]

# Add names
names(annoCharts) <- names(markers.SubclusterTypes.filtered)[setdiff(1:length(markers.SubclusterTypes.filtered), grep("Too few cells for comparison of", markers.SubclusterTypes.filtered))]

# Still some have NULL values, because no sig KEGG_PATHWAYS, remove them

annoCharts <- annoCharts[!sapply(annoCharts, is.null)]
annoCharts <- annoCharts[c(1:6,8:172)]