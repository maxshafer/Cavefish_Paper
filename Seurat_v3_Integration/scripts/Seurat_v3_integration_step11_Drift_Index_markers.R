library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

# Load and generate subsets

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

DefaultAssay(hypo.integrated) <- "RNA"

Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

Idents(hypo.integrated) <- "species"
hypo.integrated.surface <- subset(hypo.integrated, idents = "astyanax_surface")
hypo.integrated.cave <- subset(hypo.integrated, idents = "astyanax_cave")

############## Subtype Markers ##############

# Find Conserved Markers

Idents(hypo.integrated) <- "integrated_Subtype"
Idents(hypo.integrated.zeb) <- "integrated_Subtype"
Idents(hypo.integrated.ast) <- "integrated_Subtype"
Idents(hypo.integrated.surface) <- "integrated_Subtype"
Idents(hypo.integrated.cave) <- "integrated_Subtype"


conserved.markers <- lapply(levels(Idents(hypo.integrated)), function(x) FindConservedMarkers(hypo.integrated, ident.1 = x, grouping.var = "species.2", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers) <- levels(Idents(hypo.integrated))

conserved.markers.ast <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers.ast) <- levels(Idents(hypo.integrated.ast))

# Find Markers for other subsets

zebrafish.markers <- lapply(levels(Idents(hypo.integrated.zeb)), function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(zebrafish.markers) <- levels(Idents(hypo.integrated.zeb))

astyanax.markers <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(astyanax.markers) <- levels(Idents(hypo.integrated.ast))
astyanax.markers <- astyanax.markers[c(names(zebrafish.markers))]

surface.markers <- lapply(levels(Idents(hypo.integrated.surface)), function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = F, max.cells.per.ident = 1000))
names(surface.markers) <- levels(Idents(hypo.integrated.surface))

cave.markers <- lapply(levels(Idents(hypo.integrated.cave)), function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = F, max.cells.per.ident = 1000))
names(cave.markers) <- levels(Idents(hypo.integrated.cave))


############## SubclusterType Markers ##############

## Find Conserved Sub Markers
## First identify those subclusters that are species specific

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

Idents(hypo.integrated) <- "integrated_SubclusterType"
Idents(hypo.integrated.zeb) <- "integrated_SubclusterType"
Idents(hypo.integrated.ast) <- "integrated_SubclusterType"

subclusters <- levels(Idents(hypo.integrated))
index <- vector()
for(i in 1:length(levels(Idents(hypo.integrated)))) {
	zcells <- nrow(hypo.integrated.zeb@meta.data[hypo.integrated.zeb@meta.data$integrated_SubclusterType == subclusters[[i]],])
	acells <- nrow(hypo.integrated.ast@meta.data[hypo.integrated.ast@meta.data$integrated_SubclusterType == subclusters[[i]],])
	if (zcells < 3 | acells < 3){
		index[i] <- FALSE
	} else { index[i] <- TRUE}
}

conserved.markers.sub <- lapply(levels(Idents(hypo.integrated))[index], function(x) FindConservedMarkers(hypo.integrated, ident.1 = x, grouping.var = "species.2", verbose = F, max.cells.per.ident = 500))
names(conserved.markers.sub) <- levels(Idents(hypo.integrated))[index]

# Find Sub Markers for each species.2

zebrafish.markers.sub <- lapply(levels(Idents(hypo.integrated.zeb))[index], function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(zebrafish.markers.sub) <- levels(Idents(hypo.integrated.zeb))[index]

astyanax.markers.sub <- lapply(levels(Idents(hypo.integrated))[index], function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(astyanax.markers.sub) <- levels(Idents(hypo.integrated))[index]

# Find Sub Markers for each species (cave/surface)

Idents(hypo.integrated.ast) <- "integrated_SubclusterType"
Idents(hypo.integrated.surface) <- "integrated_SubclusterType"
Idents(hypo.integrated.cave) <- "integrated_SubclusterType"

subclusters <- levels(Idents(hypo.integrated.ast))
index <- vector()
for(i in 1:length(levels(Idents(hypo.integrated.ast)))) {
	scells <- nrow(hypo.integrated.surface@meta.data[hypo.integrated.surface@meta.data$integrated_SubclusterType == subclusters[[i]],])
	ccells <- nrow(hypo.integrated.cave@meta.data[hypo.integrated.cave@meta.data$integrated_SubclusterType == subclusters[[i]],])
	if (scells < 3 | ccells < 3){
		index[i] <- FALSE
	} else { index[i] <- TRUE}
}

gene.lists <- readRDS("drift_gene_lists.rds")

# conserved.markers.ast.sub <- lapply(levels(Idents(hypo.integrated.ast))[index], function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 500))
# names(conserved.markers.ast.sub) <- levels(Idents(hypo.integrated.ast))[index]

# gene.lists[[10]] <- conserved.markers.ast.sub
# saveRDS(gene.lists, file = "drift_gene_lists.rds")
# print("done conserved markers")

surface.markers.sub <- lapply(levels(Idents(hypo.integrated.ast))[index], function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(surface.markers.sub) <- levels(Idents(hypo.integrated.ast))[index]

gene.lists[[11]] <- surface.markers.sub
saveRDS(gene.lists, file = "drift_gene_lists.rds")
print("done surface markers")


cave.markers.sub <- lapply(levels(Idents(hypo.integrated.ast))[index], function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = F, max.cells.per.ident = 500))
names(cave.markers.sub) <- levels(Idents(hypo.integrated.ast))[index]

gene.lists[[12]] <- cave.markers.sub
saveRDS(gene.lists, file = "drift_gene_lists.rds")
print("done cave markers")


## Need to remove GABA_5 markers, and rename GABA_6 to GABA_5
## First remove GABA_5

for (i in 1:length(gene.lists)) {
	gene.lists[[i]] <- gene.lists[[i]][!grepl("GABA_5", names(gene.lists[[i]]))]
}

for (i in 1:length(gene.lists)) {
	names(gene.lists[[i]]) <- str_replace(names(gene.lists[[i]]), "GABA_6", "GABA_5")
}

saveRDS(gene.lists, file = "drift_gene_lists.rds")

# gene.lists <- list(conserved.markers, zebrafish.markers, astyanax.markers, conserved.markers.sub, zebrafish.markers.sub, astyanax.markers.sub, conserved.markers.ast, surface.markers, cave.markers, conserved.markers.ast.sub, surface.markers.sub, cave.markers.sub)


## Plot Marker genes

hypo.integrated@meta.data$integrated_Subtype <- factor(hypo.integrated@meta.data$integrated_Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

genes.to.plot <- lapply(gene.lists.pos[[1]], function(x) row.names(x)[1:5])
genes.to.plot <- genes.to.plot[c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia")]

DotPlot(hypo.integrated, features = unique(unlist(genes.to.plot)), group.by = "integrated_Subtype", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1), axis.title = element_blank())




