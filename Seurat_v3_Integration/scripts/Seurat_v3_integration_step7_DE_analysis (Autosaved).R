library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")
load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")


# Add species.2.Subtype and subclustertype to meta.data
hypo.integrated@meta.data$species_Subtype <- paste(hypo.integrated@meta.data$species.2, hypo.integrated@meta.data$integrated_Subtype, sep = "_")

hypo.integrated@meta.data$species_SubclusterType <- paste(hypo.integrated@meta.data$species.2, hypo.integrated@meta.data$integrated_SubclusterType, sep = "_")

# Find DE genes across species (all cells) to filter for cell type specific ones

Idents(hypo.integrated) <- "species.2"
species.markers <- FindMarkers(hypo.integrated, ident.1 = "zebrafish", ident.2 = "astyanax", only.pos = FALSE, logfc.threshold = 0.5, max.cells.per.ident = 20000)

# Find DE genes between species for all clusters (Sub and subcluster)

Idents(hypo.integrated) <- "integrated_Subtype"
Subtypes <- levels(Idents(hypo.integrated))
Idents(hypo.integrated) <- "integrated_SubclusterType"

species_Subtype <- unique(hypo.integrated@meta.data$species_Subtype)
species_SubclusterType <- unique(hypo.integrated@meta.data$species_SubclusterType)

# Find Subtype DE genes

Idents(hypo.integrated) <- "species_Subtype"
markers.Subtypes <- list()
for (i in 1:length(Subtypes)) {
	markers.Subtypes[[i]] <- FindMarkers(hypo.integrated, ident.1 = species_Subtype[grep(Subtypes[[i]], species_Subtype)][1], ident.2 = species_Subtype[grep(Subtypes[[i]], species_Subtype)][2], only.pos = FALSE, logfc.threshold = 0.5)
}

names(markers.Subtypes) <- Subtypes
# Filter by species markers

markers.Subtypes.filtered <- lapply(markers.Subtypes, function(x) x[setdiff(row.names(x), row.names(species.markers)),])

# Plots some dots

DotPlot(hypo.integrated, features = row.names(markers.Subtypes.filtered[[2]])[101:200], group.by = "integrated_Subtype", split.by = "species.2", cols = c("blue", "red")) + RotatedAxis()


# Now do the same for SubclusterTypes (more complicated because many are specific, or have too few for comparison)
SubclusterTypes <- unique(hypo.integrated@meta.data$integrated_SubclusterType)
Idents(hypo.integrated) <- "species_SubclusterType"
markers.SubclusterTypes <- list()
for (i in 1:length(SubclusterTypes)) {
	if (length(species_SubclusterType[grep(SubclusterTypes[[i]], species_SubclusterType)]) < 2){
		markers.SubclusterTypes[[i]] <- paste("Too few cells for comparison of ", SubclusterTypes[[i]], sep = "")
		print(paste("Too few cells for comparison of ", SubclusterTypes[[i]], sep = ""))
	} else {
		if (length(WhichCells(hypo.integrated, idents = species_SubclusterType[grep(SubclusterTypes[[i]], species_SubclusterType)][1])) <= 3 | length(WhichCells(hypo.integrated, idents = species_SubclusterType[grep(SubclusterTypes[[i]], species_SubclusterType)][2])) <= 3) {
		markers.SubclusterTypes[[i]] <- paste("Too few cells for comparison of ", SubclusterTypes[[i]], sep = "")
		print(paste("Too few cells for comparison of ", SubclusterTypes[[i]], sep = ""))
		} else {
			print(as.character(SubclusterTypes[[i]]))
			markers.SubclusterTypes[[i]] <- FindMarkers(hypo.integrated, ident.1 = species_SubclusterType[grep(SubclusterTypes[[i]], species_SubclusterType)][1], ident.2 = species_SubclusterType[grep(SubclusterTypes[[i]], species_SubclusterType)][2], only.pos = FALSE, logfc.threshold = 0.5)
		}
	}
}

names(markers.SubclusterTypes) <- SubclusterTypes

# Filter by species markers

markers.SubclusterTypes.filtered <- list()
for (i in 1:length(markers.SubclusterTypes)) {
	if (is.character(markers.SubclusterTypes[[i]])) {
		markers.SubclusterTypes.filtered[[i]] <- markers.SubclusterTypes[[i]]
	} else {
		markers.SubclusterTypes.filtered[[i]] <- markers.SubclusterTypes[[i]][setdiff(row.names(markers.SubclusterTypes[[i]]), row.names(species.markers)),]
	}
}

names(markers.SubclusterTypes.filtered) <- SubclusterTypes


Idents(hypo.integrated) <- "integrated_Subtype"
DotPlot(SubsetData(hypo.integrated, ident.use = "Ependymal"), features = row.names(markers.SubclusterTypes.filtered[[145]])[1:100], group.by = "integrated_SubclusterType", split.by = "species.2", cols = c("blue", "red")) + RotatedAxis()


# Save lists as rds objects

saveRDS(markers.Subtypes.filtered, file = "Zeb_Ast_markers.Subtypes.filtered.rds")
saveRDS(markers.SubclusterTypes.filtered, file = "Zeb_Ast_markers.SubclusterTypes.filtered.rds")




