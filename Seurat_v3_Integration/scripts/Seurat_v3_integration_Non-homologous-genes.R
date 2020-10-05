library(ggplot2)
library(reshape2)
library(Seurat)
library(ggpubr)
library(patchwork)
library(viridis)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

# Load normed expression datasets, which has for each integrated_SubclusterType, and for each species (zeb, ast, surface, cave)
# Load marker genes (conserved and species specific)

normed.expression <- readRDS("Normed_expression_data.rds")
markers.2 <- readRDS("drift_gene_lists_pos.rds")

hypo <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v3.rds")

Idents(hypo) <- "species.2"
hypo.integrated.zeb <- subset(hypo, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo, idents = "astyanax")

## Function

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
  data$neuronal <- ifelse(data$Subtype %in% as.character(unique(object@meta.data$integrated_Subtype)[grepl("GABA|Glut|Prdx", unique(object@meta.data$integrated_Subtype))]), "neuronal", "non-neuronal")
  data$Subclusters <- ifelse(data$Var1 %in% specific.names, "Species specific", "Conserved")
  return(data)
}

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?

prop.table <- table(hypo@meta.data$integrated_SubclusterType, hypo@meta.data$species.2)
prop.table <- prop.table[prop.table[,1] > 0 | prop.table[,2] > 0,] # I don't know why, but this still includes the removed subclusters (weird erythrocyte ones) - they are still factor levels, but also still unique entries in the meta.data??
zeb.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
ast.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

################################ Species specific gene lists (non-homologous) ################################
######### Non-homologous Genes #########

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_orthologs_corrected.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_orthologs.txt", head = TRUE)
names(mart) <- c("danio", "astyanax")

hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")

## Extract which genes don't have a homolog in the other species
ast.genes <- mart[[2]][mart[[2]]$Zebrafish.gene.stable.ID == "",]
ast.genes$Gene.name[ast.genes$Gene.name == ""] <- ast.genes$Gene.stable.ID[ast.genes$Gene.name == ""]
ast.genes <- ast.genes$Gene.name

zeb.genes0 <- mart[[1]][mart[[1]]$Cave.fish.gene.stable.ID == "",]
zeb.genes1 <- zeb.genes0
zeb.genes1$Gene.name[!(is.na(zeb.genes1$GeneID))] <- zeb.genes1$GeneID[!(is.na(zeb.genes1$GeneID))]
zeb.genes <- zeb.genes1$Gene.name

zeb.genes <- zeb.genes[zeb.genes %in% row.names(GetAssayData(hypo.zeb))]
ast.genes <- ast.genes[ast.genes %in% row.names(GetAssayData(hypo.ast))]

################################ Data, either Marker or Expressed genes ################################
######### Marker Genes #########

## Or just marker genes (from drift, plus missing subclusters for each)

zeb.markers <- markers.2[["zebrafish.markers.sub"]]
ast.markers <- markers.2[["astyanax.markers.sub"]]
zebrafish.markers.sub <- markers.2[["zeb.name.markers"]]
cavefish.markers.sub <- markers.2[["ast.name.markers"]]

zeb.markers <- c(zeb.markers[!(names(zeb.markers) %in% names(zebrafish.markers.sub))], zebrafish.markers.sub)
zeb.markers <- lapply(zeb.markers, function(x) x[x$p_val_adj < 0.0001,])
ast.markers <- c(ast.markers[!(names(ast.markers) %in% names(cavefish.markers.sub))], cavefish.markers.sub)
ast.markers <- lapply(ast.markers, function(x) x[x$p_val_adj < 0.0001,])

zeb.markers <- lapply(zeb.markers, function(x) row.names(x))
ast.markers <- lapply(ast.markers, function(x) row.names(x))

######### Expressed Genes #########

cutoff <- 10
ast.data <- normed.expression[["integrated_SubclusterType"]][["hypo.ast"]]
zeb.data <- normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]]
ast.data <- lapply(ast.data, function(x) row.names(x)[x > cutoff])
zeb.data <- lapply(zeb.data, function(x) row.names(x)[x > cutoff])

################################################# Make Plots #################################################

danio <- list()
astyanax <- list()

########### Expressed genes * Non-homologous genes

danio[[1]] <- prepData(test.data = zeb.data, test.genes = zeb.genes, object = hypo.integrated.zeb, specific.names = zeb.names)
danio[[1]]$comparison <- "Expressed Genes"
danio[[1]]$species <- "D. rerio"

astyanax[[1]] <- prepData(test.data = ast.data, test.genes = ast.genes, object = hypo.integrated.ast, specific.names = ast.names)
astyanax[[1]]$comparison <- "Expressed Genes"
astyanax[[1]]$species <- "A. mexicanus"

########### Marker genes * Non-homologous genes

danio[[2]] <- prepData(test.data = zeb.markers, test.genes = zeb.genes, object = hypo.integrated.zeb, specific.names = zeb.names)
danio[[2]]$comparison <- "Marker Genes"
danio[[2]]$species <- "D. rerio"

astyanax[[2]] <- prepData(test.data = ast.markers, test.genes = ast.genes, object = hypo.integrated.ast, specific.names = ast.names)
astyanax[[2]]$comparison <- "Marker Genes"
astyanax[[2]]$species <- "A. mexicanus"

## Save plots for figure

plot.data <- do.call(rbind, list(danio[[1]], danio[[2]], astyanax[[1]], astyanax[[2]]))
plot.data$species <- factor(plot.data$species, levels = c("D. rerio", "A. mexicanus"))
plot.data$comparison <- factor(plot.data$comparison, levels = c("Marker Genes", "Expressed Genes"))

p1 <- ggplot(plot.data[plot.data$comparison == "Marker Genes",], aes(x = Subclusters, y = value, color = Subclusters)) + geom_jitter(size = .5) + scale_color_manual(values = c("grey65", "red")) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")
p2 <- ggplot(plot.data[plot.data$comparison == "Expressed Genes",], aes(x = Subclusters, y = value, color = Subclusters)) + geom_jitter(size = .5) + scale_color_manual(values = c("grey65", "red")) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")

p1 + p2 # expressed genes with cutoff 10

## Separate by neuronal/non-neuronal

p3 <- ggplot(plot.data[plot.data$comparison == "Marker Genes",], aes(x = neuronal, y = value, color = neuronal)) + geom_jitter(size = .5) + scale_color_manual(values = c("grey65", "blue")) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")
p4 <- ggplot(plot.data[plot.data$comparison == "Expressed Genes",], aes(x = neuronal, y = value, color = neuronal)) + geom_jitter(size = .5) + scale_color_manual(values = c("grey65", "blue")) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")

p3 + p4 # Expressed genes has clear effect (zeb neuronal, ast non-neuronal)

# Plot by Subtype

p5 <- ggplot(plot.data[plot.data$comparison == "Marker Genes",], aes(x = Subtype, y = value, color = Subtype)) + geom_jitter(size = .5) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")
p6 <- ggplot(plot.data[plot.data$comparison == "Expressed Genes",], aes(x = Subtype, y = value, color = Subtype)) + geom_jitter(size = .5) + geom_boxplot(color = "black", outlier.shape = NA, fill = NA) + theme_classic() + theme(axis.text.x = element_text(angle = 45, size = 8, hjust = 1, vjust = 1), axis.text = element_text(size = 8), axis.title = element_text(size = 10)) + stat_compare_means(size = 3) + guides(color = F) + ylab("% non-homologous genes") + xlab("") + facet_wrap(~comparison + species, scales = "free_y")

p5 / p6 # marker genes looks better

nhgenes <- p1 + p2 + p3 + p4

nhgenes <- nhgenes + plot_layout(width = unit(c(40), c("mm")), height = unit(c(40), c("mm")))
dev.new()
nhgenes









## Plot sox1a featureplots for figure
Idents(hypo) <- "integrated_Subtype"
hypo.neuronal <- subset(hypo, idents = as.character(unique(hypo@meta.data$integrated_Subtype)[grepl("GABA|Glut|Prdx", unique(hypo@meta.data$integrated_Subtype))]))



embed <- Embeddings(hypo.neuronal, reduction = "tsne")
meta <- hypo.neuronal@meta.data[,c("species", "species.2", "integrated_Subtype", "integrated_SubclusterType")]
meta <- meta[row.names(meta) %in% row.names(embed),]

tsne.data <- cbind(embed, meta)
tsne.data$cell_id <- row.names(tsne.data)
tsne.data$highlight <- "absent"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("Glut_7_2")] <- "Glut_7_2 - vipb+"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("GABA_0_7")] <- "vip+"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("GABA_2_1")] <- "vip+"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("GABA_2_3")] <- "vip+"
tsne.data$highlight[tsne.data$integrated_SubclusterType == c("GABA_2_8")] <- "vip+"
tsne.data$vip <- GetAssayData(hypo.neuronal)["vip",]
tsne.data$vipb <- GetAssayData(hypo.neuronal)["vipb",]
tsne.data <- tsne.data[order(tsne.data$vipb, decreasing = F)[order(tsne.data$highlight, decreasing = F)],]



big.plot <- ggplot(tsne.data, aes(x = basetsne_1, y = basetsne_2, fill = vip, color = highlight)) + theme_classic() + geom_point(size = 1, shape = 21) + scale_fill_gradient(high = "blue", low = "grey90") + scale_color_manual(values = c("transparent", "yellow", "red")) + theme(axis.title = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), axis.text = element_blank()) + geom_rect(xmin = 4, xmax = 8, ymin = 22, ymax = 26, fill = "transparent", color = "black")

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

