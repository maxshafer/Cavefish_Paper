library(Seurat)
library(tidyr)
library(dendextend)
library(dplyr)
library(ape)
library(circlize)
library(ggplot2)
library(scales)
library(tibble)

# Load seurat object (astmex)
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/")
load("AstMex_64k.Robj")


# Find Subtype species markers

hypo <- SetAllIdent(hypo, id = "Subtype")

subtypes <- levels(hypo@ident)

hypo <- SetAllIdent(hypo, id = "species_subtype")

markers <- list()i
for (i in 1:length(subtypes)){
	markers[[i]] <- FindMarkers(hypo, ident.1 = paste("astyanax_surface_", subtypes[[i]], sep = ""), ident.2 = paste("astyanax_cave_", subtypes[[i]], sep = ""), only.pos = FALSE, logfc.threshold = 0.5)
}

markers <- lapply(subtypes, function(x) FindMarkers(hypo, ident.1 = paste("astyanax_surface_", x, sep = ""), ident.2 = paste("astyanax_cave_", x, sep = ""), only.pos = FALSE, logfc.threshold = 0.5))

names(markers) <- subtypes
saveRDS(markers, file = "AstMex_Hypo_markers.species.Subtype.list.rds")



# Find SubclusterType species markers

hypo <- SetAllIdent(hypo, id = "SubclusterType")

subclustertypes <- levels(hypo@ident)

hypo <- SetAllIdent(hypo, id = "species_subcluster")


## If they exist
## If they have more then x number per cluster
## Calculate even for super uneven subclusters, decided to not show them on graph later

markers <- list()
for (i in 51:length(subclustertypes)){
	if(paste("astyanax_cave_", subclustertypes[[i]], sep = "") %in% levels(hypo@ident) & paste("astyanax_surface_", subclustertypes[[i]], sep = "") %in% levels(hypo@ident)) {
		if(length(WhichCells(hypo, ident = paste("astyanax_cave_", subclustertypes[[i]], sep = ""))) > 3 & length(WhichCells(hypo, ident = paste("astyanax_surface_", subclustertypes[[i]], sep = ""))) > 3) {
			markers[[i]] <- FindMarkers(hypo, ident.1 = paste("astyanax_surface_", subclustertypes[[i]], sep = ""), ident.2 = paste("astyanax_cave_", subclustertypes[[i]], sep = ""), only.pos = FALSE, logfc.threshold = 0.5)
		}
	} else {
		markers[[i]] <- "morph specific cell type"
	}
}


names(markers) <- subclustertypes
saveRDS(markers, file = "AstMex_Hypo_markers.species.SubclusterType.list.rds")