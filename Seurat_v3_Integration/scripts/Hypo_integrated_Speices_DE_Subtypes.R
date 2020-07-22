library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(tidyr)
library(dendextend)
library(dplyr)
library(ape)
library(circlize)
library(ggplot2)
library(scales)
library(tibble)

# Load seurat object (integrated)
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")
hypo <- hypo.integrated

# Find Subtype species markers

Idents(hypo) <- "integrated_Subtype"
subtypes <- levels(Idents(hypo))

Idents(hypo) <- "integrated_Subtype_species"

markers <- list()
for (i in 1:length(subtypes)){
	markers[[i]] <- FindMarkers(hypo, ident.1 = paste("zebrafish_", subtypes[[i]], sep = ""), ident.2 = paste("astyanax_", subtypes[[i]], sep = ""), only.pos = FALSE, logfc.threshold = 0.5)
}

markers <- lapply(subtypes, function(x) FindMarkers(hypo, ident.1 = paste("astyanax_surface_", x, sep = ""), ident.2 = paste("astyanax_cave_", x, sep = ""), only.pos = FALSE, logfc.threshold = 0.5))

names(markers) <- subtypes
saveRDS(markers, file = "Hypo_integrated_markers.species.Subtype.list.rds")



# Find SubclusterType species markers

Idents(hypo) <- "integrated_SubclusterType"
subclustertypes <- levels(Idents(hypo))

Idents(hypo) <- "integrated_SubclusterType_species"


## If they exist
## If they have more then x number per cluster
## Calculate even for super uneven subclusters, decided to not show them on graph later

markers <- list()
for (i in 1:50){
	if(paste("zebrafish_", subclustertypes[[i]], sep = "") %in% levels(Idents(hypo)) & paste("astyanax_", subclustertypes[[i]], sep = "") %in% levels(Idents(hypo))) {
		if(length(WhichCells(hypo, ident = paste("zebrafish_", subclustertypes[[i]], sep = ""))) > 3 & length(WhichCells(hypo, ident = paste("astyanax_", subclustertypes[[i]], sep = ""))) > 3) {
			markers[[i]] <- FindMarkers(hypo, ident.1 = paste("zebrafish_", subclustertypes[[i]], sep = ""), ident.2 = paste("astyanax_", subclustertypes[[i]], sep = ""), only.pos = FALSE, logfc.threshold = 0.5)
		}
	} else {
		markers[[i]] <- "species specific cell type"
	}
}


names(markers) <- subclustertypes
saveRDS(markers, file = "Hypo_integrated_markers.species.SubclusterType.list.rds")