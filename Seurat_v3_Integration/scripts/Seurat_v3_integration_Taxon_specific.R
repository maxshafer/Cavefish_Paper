library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(Seurat)
library(ggpubr)
library(patchwork)
library(SuperExactTest)


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")
hypo <- hypo.integrated

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?

prop.table <- table(hypo@meta.data$integrated_SubclusterType, hypo@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$Zeb_specific <- ifelse(prop.table$zebrafish > .9, "yes", "no")
prop.table$Ast_specific <- ifelse(prop.table$zebrafish < .1, "yes", "no")

zeb.names <- row.names(prop.table[prop.table$Zeb_specific == "yes",])
ast.names <- row.names(prop.table[prop.table$Ast_specific == "yes",])

################################ Species specific gene lists (taxon, non-homologous, jialin's) ################################

######### Taxon-specific duplications ###########

## Load taxon specific gene lists
### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart.zeb <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Otophysi_Paralogs/mart_export_Danio.txt", sep = "\t", head = TRUE)
mart.ast <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Otophysi_Paralogs/mart_export_Astyanax.txt", sep = "\t", head = TRUE)
# Node levels for phylo tree for the origin of duplicated gene pairs
nodes.zeb <- levels(mart.zeb$Paralogue.last.common.ancestor.with.Zebrafish)
nodes.ast <- levels(mart.ast$Paralogue.last.common.ancestor.with.Mexican.tetra)

# Subset mart for each of the node levels and extract unique gene.names
mart.zeb <- lapply(nodes.zeb, function(x) as.character(unique(mart.zeb[mart.zeb$Paralogue.last.common.ancestor.with.Zebrafish == x, 2])))
names(mart.zeb) <- nodes.zeb
mart.ast <- lapply(nodes.ast, function(x) as.character(unique(mart.ast[mart.ast$Paralogue.last.common.ancestor.with.Mexican.tetra == x, 2])))
names(mart.ast) <- nodes.ast

######### Non-homologous Genes #########

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_danio.txt", head = TRUE)
mart[[2]] <- read.csv("~/Downloads/mart_export_AstMex102.txt", stringsAsFactors = F)
names(mart) <- c("danio", "astyanax")

## Extract which genes don't have a homolog in the other species
ast.genes <- mart[[2]][mart[[2]]$ZFIN.ID == "",]
ast.genes$Gene.name[ast.genes$Gene.name == ""] <- ast.genes$Gene.stable.ID[ast.genes$Gene.name == ""]
ast.genes <- ast.genes$Gene.name

zeb.genes <- mart[[1]]$Gene.name[mart[[1]]$Mexican.tetra.gene.stable.ID == ""]

load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj")
load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_65k.Robj")

zeb.genes <- unique(zeb.genes[zeb.genes %in% row.names(hypo.zeb@data)])
ast.genes <- ast.genes[ast.genes %in% row.names(hypo.ast@data)]

######### Jialin Genes #########

Idents(hypo.integrated) <- "species.2"
DefaultAssay(hypo.integrated) <- "RNA"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

## Zeb specific genes, from Jialin
library(org.Dr.eg.db)

jialin.zeb <- read.csv("ageZf.txt", header = F, sep = "\t")
jialin.zeb$V3 <- mapIds(org.Dr.eg.db, keys = as.character(jialin.zeb$V1), keytype = "ENSEMBL", column = "SYMBOL")
jialin.zeb <- jialin.zeb[jialin.zeb$V3 %in% row.names((GetAssayData(hypo.integrated.zeb, slot = "counts"))),] # This doesn't affect the % in the data, but % of the data in this
names <- unique(jialin.zeb$V2)
jialin.zeb <- lapply(names, function(x) jialin.zeb[jialin.zeb$V2 == x, 3])
names(jialin.zeb) <- names

################################ Data, either Marker or Expressed genes ################################

######### Marker Genes #########

## Or just marker genes (from drift, plus missing subclusters for each)

markers.2 <- readRDS("drift_gene_lists.rds")

Idents(hypo.integrated.zeb) <- "integrated_SubclusterType"
Idents(hypo.integrated.ast) <- "integrated_SubclusterType"

zebrafish.markers.sub <- lapply(zeb.names, function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(zebrafish.markers.sub) <- zeb.names

cavefish.markers.sub <- lapply(ast.names, function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(cavefish.markers.sub) <- ast.names

zeb.markers <- markers.2[[5]]
ast.markers <- markers.2[[6]]

zeb.markers <- c(zeb.markers[!(names(zeb.markers) %in% names(zebrafish.markers.sub))], zebrafish.markers.sub)
zeb.markers <- lapply(zeb.markers, function(x) x[x$avg_logFC > 0.5,])
ast.markers <- c(ast.markers[!(names(ast.markers) %in% names(cavefish.markers.sub))], cavefish.markers.sub)
ast.markers <- lapply(ast.markers, function(x) x[x$avg_logFC > 0.5,])

zeb.markers <- lapply(zeb.markers, function(x) row.names(x))
ast.markers <- lapply(ast.markers, function(x) row.names(x))

######### Expressed Genes #########

# Load normed expression datasets, which has for each integrated_SubclusterType, and for each species (zeb, ast, surface, cave)
normed.expression <- readRDS("Normed_expression_data.rds")

str(normed.expression, max.level = 2)

# ## Remove and resave normed.expression
# ## remove GABA_5, and rename GABA_6 clusters for integrated Sutype and SubclusterType

# ast.data <- normed.expression[["integrated_Subtype"]][["hypo.ast"]][!grepl("GABA_5", names(normed.expression[["integrated_Subtype"]][["hypo.ast"]]))]
# zeb.data <- normed.expression[["integrated_Subtype"]][["hypo.zeb"]][!grepl("GABA_5", names(normed.expression[["integrated_Subtype"]][["hypo.zeb"]]))]

# names(zeb.data)[grep("GABA_6", names(zeb.data))] <- c("GABA_5")
# names(ast.data)[grep("GABA_6", names(ast.data))] <- c("GABA_5")

# normed.expression[["integrated_Subtype"]][["hypo.ast"]] <- ast.data
# normed.expression[["integrated_Subtype"]][["hypo.zeb"] <- zeb.data

# ast.data <- normed.expression[["integrated_SubclusterType"]][["hypo.ast"]][!grepl("GABA_5", names(normed.expression[["integrated_SubclusterType"]][["hypo.ast"]]))]
# zeb.data <- normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]][!grepl("GABA_5", names(normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]]))]

# names(zeb.data)[grep("GABA_6", names(zeb.data))] <- c("GABA_5_0", "GABA_5_1", "GABA_5_2")
# names(ast.data)[grep("GABA_6", names(ast.data))] <- c("GABA_5_0", "GABA_5_1", "GABA_5_2")

# normed.expression[["integrated_SubclusterType"]][["hypo.ast"]] <- ast.data
# normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]] <- zeb.data

# saveRDS(normed.expression , "Normed_expression_data.rds")

# What I want is normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]][zeb.names] and normed.expression[["integrated_SubclusterType"]][["hypo.ast"]][ast.names]
# lists of averaged expression for each gene in each cell type



cutoff <- 2
ast.data <- normed.expression[["integrated_SubclusterType"]][["hypo.ast"]]
zeb.data <- normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]]
ast.data <- lapply(ast.data, function(x) row.names(x)[x > cutoff])
zeb.data <- lapply(zeb.data, function(x) row.names(x)[x > cutoff])


################################################# Make Plots #################################################

prepData <- function(test.data = test.data, test.genes = test.genes, object = object, specific.names = specific.names) {
	if("Bilateria" %in% names(test.genes)) {
		data <- lapply(test.data, function(x) lapply(test.genes, function(y) (length(x[x %in% y])/length(x))*100))
		data <- lapply(data, function(x) unlist(x))
		data <- do.call(rbind, data)
		colnames(data)[1] <- "All"
		data <- as.data.frame(melt(data))
		data$Var2 <- factor(data$Var2, levels = c("All", "Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Otophysa", "Characiphysae", "Characoidei", "Astyanax mexicanus", "Danio rerio"))
	} else {
		data <- lapply(test.data, function(x) length(x[x %in% test.genes])/length(x)*100)
		data <- as.data.frame(unlist(data))
		colnames(data)[1] <- "value"
		data$Var1 <- row.names(data)
		}
	data$Subtype <- object@meta.data$integrated_Subtype[match(data$Var1, object@meta.data$integrated_SubclusterType)]
	data$Neuronal <- object@meta.data$neuronal[match(data$Var1, object@meta.data$integrated_SubclusterType)]
	data$Subclusters <- ifelse(data$Var1 %in% specific.names, "Species specific", "Conserved")
	return(data)
}

danio <- list()
astyanax <- list()

########### Expressed genes * Non-homologous genes

danio[[1]] <- prepData(test.data = zeb.data, test.genes = zeb.genes, object = hypo.integrated.zeb, specific.names = zeb.names)
astyanax[[1]] <- prepData(test.data = ast.data, test.genes = ast.genes, object = hypo.integrated.ast, specific.names = ast.names)

########### Marker genes * Non-homologous genes

danio[[2]] <- prepData(test.data = zeb.markers, test.genes = zeb.genes, object = hypo.integrated.zeb, specific.names = zeb.names)
astyanax[[2]] <- prepData(test.data = ast.markers, test.genes = ast.genes, object = hypo.integrated.ast, specific.names = ast.names)

########### Expressed genes * Taxon specific genes

danio[[3]] <- prepData(test.data = zeb.data, test.genes = mart.zeb, object = hypo.integrated.zeb, specific.names = zeb.names)
astyanax[[3]] <- prepData(test.data = ast.data, test.genes = mart.ast, object = hypo.integrated.ast, specific.names = ast.names)

########### Marker genes * Taxon specific genes

danio[[4]] <- prepData(test.data = zeb.markers, test.genes = mart.zeb, object = hypo.integrated.zeb, specific.names = zeb.names)
astyanax[[4]] <- prepData(test.data = ast.markers, test.genes = mart.ast, object = hypo.integrated.ast, specific.names = ast.names)

########### Marker genes and Expressed genes * Jialin specific genes

danio[[5]] <- prepData(test.data = zeb.data, test.genes = jialin.zeb, object = hypo.integrated.zeb, specific.names = zeb.names)
danio[[6]] <- prepData(test.data = zeb.markers, test.genes = jialin.zeb, object = hypo.integrated.zeb, specific.names = zeb.names)



## Save plots for figure

danio[[1]]$comparison <- "Marker Genes"
danio[[2]]$comparison <- "Expressed Genes"
danio[[1]]$species <- "D. rerio"
danio[[2]]$species <- "D. rerio"

astyanax[[1]]$comparison <- "Marker Genes"
astyanax[[2]]$comparison <- "Expressed Genes"
astyanax[[1]]$species <- "A. mexicanus"
astyanax[[2]]$species <- "A. mexicanus"

plot.data <- do.call(rbind, list(danio[[1]], danio[[2]], astyanax[[1]], astyanax[[2]]))
plot.data$species <- factor(plot.data$species, levels = c("D. rerio", "A. mexicanus"))
plot.data$comparison <- factor(plot.data$comparison, levels = c("Marker Genes", "Expressed Genes"))
plot.data$Subtype <- factor(plot.data$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

p1 <- ggplot(plot.data[plot.data$comparison == "Marker Genes",], aes(x = Subclusters, y = value, color = Subclusters)) + geom_jitter(size = .5) + scale_color_manual(values = c("grey65", "red")) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")


p2 <- ggplot(plot.data[plot.data$comparison == "Marker Genes",], aes(x = Neuronal, y = value, color = Neuronal)) + geom_jitter(size = .5) + scale_color_manual(values = c("grey65", "darkblue")) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")

p3 <- ggplot(plot.data[plot.data$comparison == "Marker Genes",], aes(x = Subtype, y = value, color = Subtype)) + geom_jitter(size = 1) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y", nrow = 2)

p4 <- p1 / p2
(p4 | p3) + plot_layout(widths = c(1,3))

## SuperExactTest

zeb.markers.2 <- lapply(zeb.markers, function(x) x[x %in% zeb.genes])
ast.markers.2 <- lapply(ast.markers, function(x) x[x %in% ast.genes])

zeb.markers.3 <- zeb.markers[zeb.names]
ast.markers.3 <- ast.markers.2[ast.names]

length.gene.sets <- sapply(zeb.markers.3,length)

total <- nrow(hypo.zeb@data)

# Just neuronal zebrafish clusters
res=supertest(zeb.markers.3[c(2,3,5:9)], n=total)
	
plot(res, Layout="landscape", degree = 2:11, sort.by="size", keep=FALSE, show.elements=TRUE, elements.cex=0.7, show.fold.enrichment = TRUE, elements.list=subset(summary(res)$Table,Observed.Overlap <= 6), show.expected.overlap=TRUE,expected.overlap.style="hatchedBox", color.expected.overlap='red')

## Plot sox1a featureplots for figure

embed <- GetCellEmbeddings(hypo.zeb, reduction.type = "tsne")
meta <- hypo.zeb@meta.data[,c(18,20,22,23)]

tsne.data <- cbind(embed, meta)
tsne.data$cell_id <- row.names(tsne.data)
tsne.data$highlight <- "absent"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("Progenitors_9")] <- "Progenitor"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("Glut_0_3")] <- "D. rerio specific subcluster"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("GABA_0_15")] <- "D. rerio specific subcluster"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("Prdx1_Positive_2")] <- "D. rerio specific subcluster"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("Prdx1_Positive_10")] <- "D. rerio specific subcluster"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("Prdx1_Positive_11")] <- "D. rerio specific subcluster"
tsne.data$sox1a <- hypo.zeb@data["sox1a",]
tsne.data <- tsne.data[order(tsne.data$sox1a, decreasing = F)[order(tsne.data$highlight, decreasing = F)],]



big.plot <- ggplot(tsne.data, aes(x = tSNE_1, y = tSNE_2, fill = sox1a, color = highlight)) + geom_point(size = 1, shape = 21) + scale_fill_gradient(high = "blue", low = "grey90") + scale_color_manual(values = c("transparent", "yellow", "red")) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) + geom_rect(xmin = 4, xmax = 8, ymin = 22, ymax = 26, fill = "transparent", color = "black")

inset <- ggplot(tsne.data, aes(x = tSNE_1, y = tSNE_2, fill = sox1a, color = highlight)) + geom_point(size = 3, shape = 21) + scale_fill_gradient(high = "blue", low = "grey90") + scale_color_manual(values = c("transparent", "yellow", "red")) + guides(color = F, fill = F) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) + xlim(c(4,8)) + ylim(c(22,26)) + ggtitle("Progenitors_9") + geom_rect(xmin = 4, xmax = 8, ymin = 22, ymax = 26, fill = "transparent", color = "black")

layout <- c(area(t = 1, b = 16, l = 2, r = 12), area(t = 1, b = 3, l = 11, r = 12))
png("sox1a_inset_test.png", width = 3, height = 4)
big.plot + inset + plot_layout(design = layout)
ggsave("sox1a_inset_test.png", width = 12, height = 9)
dev.off()


hypo.zeb <- SetAllIdent(hypo.zeb, id = "integrated_Subtype")

dplot <- DotPlot(SubsetData(hypo.zeb, ident.use = "Progenitors"), genes.plot = c("sox1a"), group.by = "integrated_SubclusterType", do.return = T, plot.legend = T) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 8))

dplot + guides(size = F)

DotPlot(SubsetData(hypo.zeb, ident.use = "Progenitors"), genes.plot = sort(zeb.markers[["Progenitors_9"]])[sort(zeb.markers[["Progenitors_9"]]) %in% zeb.genes], group.by = "SubclusterType", do.return = T) + coord_flip() + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 8))


FeaturePlot(hypo.zeb, features = c("sox1a"))











hypo.zeb.2 <- UpdateSeuratObject(hypo.zeb)

FeaturePlot(hypo.zeb.2, features = c("gpsm2"), order = T) + NoAxes() +NoLegend()

Idents(hypo.zeb.2) <- "integrated_SubclusterType"

DimPlot(hypo.zeb.2, reduction.use = "tsne", group.by = "integrated_SubclusterType", cells.highlight = WhichCells(hypo.zeb.2, idents = c("Progenitors_9", "Glut_0_3", "GABA_0_15", "Prdx1_Positive_2", "Prdx1_Positive_10", "Prdx1_Positive_11")), no.legend = T, do.label = T) + NoAxes() + NoLegend()


DotPlot(hypo.zeb, genes.plot = zeb.markers.2[["Progenitors_9"]], group.by = "integrated_SubclusterType", do.return = T) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + coord_flip()

