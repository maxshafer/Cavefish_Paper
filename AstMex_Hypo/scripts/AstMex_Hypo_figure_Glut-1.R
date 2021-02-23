library(Seurat)
library(stringr)
library(dplyr)
library(patchwork)
library(viridis)
library(RDAVIDWebService)
library(ggtree)
library(ape)
library(phylogram)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")
hypo.ast <- readRDS("AstMex_63k_vR.rds")

Idents(hypo.ast) <- "Cluster"

# Find morph specific clusters
prop.table <- table(hypo@meta.data$Subcluster, hypo@meta.data$species)
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Neuronal_03 has two cave-specific subclusters (stress response?)

neuronal_03 <- subset(hypo.ast, idents = "Neuronal_03")

neuronal_03 <- FindVariableFeatures(neuronal_03, selection.method = "mvp")
neuronal_03 <- ScaleData(object = neuronal_03, features = VariableFeatures(neuronal_03))
neuronal_03 <- RunPCA(object = neuronal_03, features = VariableFeatures(neuronal_03), npcs = 100, set.seed = 0)
# ElbowPlot(object = neuronal_03, ndims = 100) # 25 PCs looks good
neuronal_03 <- RunTSNE(object = neuronal_03, reduction = "pca", dims = 1:15, tsne.method = "Rtsne", reduction.name = "tsne", reduction.key = "tsne_", seed.use = 1, check_duplicates = F)

Idents(neuronal_03) <- "Subcluster"
neuronal_03 <- BuildClusterTree(neuronal_03, features = VariableFeatures(neuronal_03))

saveRDS(neuronal_03, file = "AstMex_neuronal_03.rds")


dendro.plot <- ggtree(as.phylo(as.dendrogram(Tool(neuronal_03, slot = "BuildClusterTree")))) + theme_tree(bgcolor = NA) #+ geom_tiplab(offset = 0)

## Make proportion plot
cols0 <- c("#FDE725FF", "#22A884FF")

prop.table <- table(neuronal_03@meta.data$Subcluster, neuronal_03@meta.data$species)[unique(neuronal_03@meta.data$Subcluster),]
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$cell_type <- row.names(prop.table)
prop.table <- reshape::melt(prop.table)
prop.table$cell_type <- factor(prop.table$cell_type, levels = as.phylo(as.dendrogram(Tool(neuronal_03, slot = "BuildClusterTree")))[[4]])


prop.plots <- ggplot(prop.table, aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity")
prop.plots <- prop.plots + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_blank()) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + ylab("Sample Cluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Origin (cave or surface)"))  + scale_fill_manual(values = cols0)

# Combine dendrogram and proportion plots
dend.plot <- dendro.plot + prop.plots + plot_layout(ncol = 2, widths = c(3,1))


## Find Markers

# cell.types <- unique(neuronal_03@meta.data$species_Subcluster)
# cell.types <- cell.types[table(Idents(neuronal_03)) > 3]
Idents(neuronal_03) <- "Subcluster"
# Idents(neuronal_03) <- factor(Idents(neuronal_03), levels = sort(unique(neuronal_03@meta.data$species_Subcluster)))
species_subtype_markers <- FindAllMarkers(neuronal_03, only.pos = T)
markers <- species_subtype_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

## Make Plots

cols0 <- c("#FDE725FF", "#22A884FF")
cols3 <- c("#771155", "#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77", "#771122", "#AA4455", "#DD7788")


neuronal_03.morph <- DimPlot(object = neuronal_03, group.by = "species", reduction = "tsne", pt.size = .25, label = FALSE, cols = cols0) + NoAxes() + theme(legend.position = c(0.8,0.9), legend.background = element_blank()) + guides(color = guide_legend(ncol = 1, override.aes = list(size = 2)))
neuronal_03.subtype <- DimPlot(object = neuronal_03, group.by = "Cluster", reduction = "tsne", pt.size = .25, label = TRUE) + NoLegend() + NoAxes()
neuronal_03.subcluster <- DimPlot(object = neuronal_03, group.by = "Subcluster", reduction = "tsne", pt.size = .25, label = TRUE) + NoLegend() + NoAxes()
neuronal_03.orig <- DimPlot(object = neuronal_03, group.by = "orig.ident", reduction = "tsne", pt.size = .25, label = FALSE, cols = cols3) + NoAxes() + theme(legend.position = c(0.8,0.9), legend.background = element_blank()) + guides(color = guide_legend(ncol = 2, override.aes = list(size = 5))) + scale_colour_manual(values = cols3)

gene <- FeaturePlot(neuronal_03, features = c("hspb1"), reduction = "tsne", pt.size = .25) + NoAxes() + ggtitle("") + theme(legend.position = c(0.8,0.9), legend.background = element_blank())

## Make dot plot
dot.plot <- DotPlot(neuronal_03, features = unique(markers$gene), group.by = "Subcluster", scale.max = 124) + NoLegend() + theme(axis.text = element_text(size = 8), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_color_viridis()


## Extract marker genes for DAVID analysis

marker.genes <- species_subtype_markers[species_subtype_markers$cluster == "Neuronal_03_4", "gene"]
marker.genes <- select(org.Dr.eg.db, marker.genes, "ENTREZID", "SYMBOL")[,2]
david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))
result <- addList(david, marker.genes, idType = "ENTREZ_GENE_ID", listName = "Neuronal_03_4", listType = "Gene")
annoCharts <- getFunctionalAnnotationChart(david)

anno.data <- annoCharts[,c(1,2,3,5,7,8,9,10,11,12,13)]

# Make DAVID plot

david.plot <- ggplot(anno.data[anno.data$PValue < 0.01,], aes(x=Term, y=Count, fill=log(Bonferroni))) + geom_bar(stat="identity") + coord_flip() + theme_classic() + theme(axis.text.y = element_blank())

david.plot <- david.plot + geom_text(aes(x=Term, y=Count+0.5, label = Term), hjust = 0) + theme(legend.position = c(0.75,0.75), legend.background = element_blank())

## Combine plots using patchwork

tsnes <- (neuronal_03.subcluster + neuronal_03.morph + plot_layout(ncol = 2)) 

row1 <- tsnes + dend.plot + plot_layout(ncol = 3, widths = c(1,1,2))

row2 <- dot.plot + plot_spacer() + plot_layout(ncol = 2, widths = c(2.75,1.25))

row3 <- david.plot + gene + plot_layout(ncol = 2, widths = c(3,1))

((row1 / row2) / row3) + plot_layout(nrow = 3)


## Small Sankey plot comes from separate script


## Make small figure for galnain and oxytocin cluster expression

galn <- subset(hypo.ast, idents = "Neuronal_07")
oxt <- subset(hypo.ast, idents = "Neuronal_19")

dot.galn <- DotPlot(galn, features = c("galn"), group.by = "morph", scale.min = 30, scale.max = 80, dot.scale = 4) 
dot.galn <- dot.galn + theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8)) + scale_color_viridis_c(limits = c(-1.5,1.5), option = "B")

dot.oxt <- DotPlot(oxt, features = c("oxt", "avp", "ENSAMXG00000021172"), group.by = "morph", scale.min = 30, scale.max = 80, dot.scale = 4) 
dot.oxt <- dot.oxt + theme(axis.title = element_blank(), axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_blank()) + scale_color_viridis_c(limits = c(-1.5,1.5), option = "B")

dev.new()
dot.galn + dot.oxt + plot_layout(guides = "collect", widths = unit(c(7,13), "mm"), height = unit(c(35,35), "mm"))




