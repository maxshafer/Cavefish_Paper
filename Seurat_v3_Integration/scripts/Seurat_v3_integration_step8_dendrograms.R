library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(dendextend)
library(ape)
library(circlize)
library(ggplot2)
library(scales)
library(phytools)
library(ggtree)
library(ggplotify)

# Load subsets
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")
load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

hypo.integrated@meta.data$species_morph[hypo.integrated@meta.data$species == "zebrafish"] <- "zebrafish"

DefaultAssay(hypo.integrated) <- "integrated"

# Make Dendrograms as a new list

dendrograms <- list()

Idents(hypo.integrated) <- "integrated_Subtype"
hypo.integrated <- BuildClusterTree(hypo.integrated)
dendrograms[[1]] <- Tool(hypo.integrated, slot = "BuildClusterTree")

Idents(hypo.integrated) <- "integrated_SubclusterType"
hypo.integrated <- BuildClusterTree(hypo.integrated)
dendrograms[[2]] <- Tool(hypo.integrated, slot = "BuildClusterTree")

Idents(hypo.integrated) <- "integrated_Subtype_species"
hypo.integrated <- BuildClusterTree(hypo.integrated)
dendrograms[[3]] <- Tool(hypo.integrated, slot = "BuildClusterTree")

Idents(hypo.integrated) <- "integrated_SubclusterType_species"
hypo.integrated <- BuildClusterTree(hypo.integrated)
dendrograms[[4]] <- Tool(hypo.integrated, slot = "BuildClusterTree")


############ Make dendrograms for each cluster and subcluster for species/morphs

# Add meta data which includes zebrafish and all 4 morphs

hypo.integrated@meta.data$species.3 <- ifelse(hypo.integrated@meta.data$species.2 %in% "zebrafish", "zebrafish", ifelse(hypo.integrated@meta.data$species %in% "astyanax_surface", "astyanax_surface", paste("astyanax",hypo.integrated@meta.data$species_morph, sep = "_")))

## Function for calculating dendrograms by species for each cluster/subcluster

denDro <- function(object = hypo.integrated, cluster = cluster, idents = "integrated_Subtype", assay = "RNA") {
	Idents(object) <- idents
	object <- subset(object, idents = cluster, features = genes)
	Idents(object) <- "species.3"
	DefaultAssay(object) <- assay
	object <- BuildClusterTree(object, features = genes)
	return(Tool(object, slot = "BuildClusterTree"))
}

## For clusters/cell types
 
dendrograms.subtype <- lapply(unique(hypo.integrated@meta.data$integrated_Subtype), function(x) denDro(object = hypo.integrated, cluster = x, idents = "integrated_Subtype", assay = "RNA"))

names(dendrograms.subtype) <- unique(hypo.integrated@meta.data$integrated_Subtype)

dendrograms[[5]] <- dendrograms.subtype


# # Plot all trees as facets
# class(dendrograms.subtype) <- "multiPhylo"
# ggtree(dendrograms.subtype) + facet_wrap(~.id, scale = "free") + geom_tiplab()

densityTree(dendrograms.subtype, type = "cladogram", use.edge.length = F, nodes = "intermediate", show.axis = F, alpha = 0.1, color = "black")

## For Subclusters

SubclusterTypes <- unique(hypo.integrated@meta.data$integrated_SubclusterType)
index <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species.3)
index <- t(apply(index, 1, function(x) x > 2))
index2 <- apply(index, 1, function(x) length(x[x]) > 4)
index <- apply(index, 1, function(x) length(x[x]) > 2)

dendrograms.subclustertype <- lapply(SubclusterTypes[SubclusterTypes %in% names(index)[index]], function(x) denDro(object = hypo.integrated, cluster = x, idents = "integrated_SubclusterType", assay = "RNA"))

names(dendrograms.subclustertype) <- SubclusterTypes[SubclusterTypes %in% names(index)[index]]

dendrograms[[6]] <- dendrograms.subclustertype

# # Plot all trees as facets
# class(dendrograms.subclustertype) <- "multiPhylo"
# ggtree(dendrograms.subclustertype) + facet_wrap(~.id, scale = "free") + geom_tiplab()

# Plot density tree for subclusters which are shared by all morphs + zebrafish
densityTree(dendrograms.subclustertype[names(dendrograms.subclustertype) %in% names(index)[index & index2]], type = "cladogram", use.edge.length = F, nodes = "intermediate", show.axis = F, alpha = 0.02, color = "black")



##### Save dendrograms

# Name and as.dendrogram the list


names(dendrograms) <- c("integrated_Subtype", "integrated_SubclusterType", "species_Subtype", "species_SubclusterType", "Subtype_species-species-morph", "SubclusterType_species-species-morph")

dendrograms <- lapply(dendrograms, function(x) as.dendrogram(x))


saveRDS(dendrograms, file = "Zeb_Ast_dendrograms_integrated.rds")

## Don't know what this is for, but I'll keep it for now
## It makes heatmaps for dendrogram distances - like the correlation plot I have in Figure 2

phylo <- list()
phylo <- lapply(dendrograms, function(x) cophenetic.phylo(as.phylo(x)))

mphylo <- melt(phylo[[3]])
mphylo <- mphylo[grep("zebrafish", mphylo$Var1),]
mphylo <- mphylo[grep("astyanax", mphylo$Var2),]
ggplot(data = mphylo, aes(x=Var1, y=Var2, fill = log(value))) + geom_tile() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

get_lower_tri <- function(mphylo) {
	mphylo[upper.tri(mphylo)] <- NA
	return(mphylo)
}

get_upper_tri <- function(mphylo) {
	mphylo[lower.tri(mphylo)] <- NA
	return(mphylo)
}

upper_tri <- get_upper_tri(phylo[[3]])

melted_phylo <- melt(upper_tri, na.rm = TRUE)


ggplot(data = melted_phylo, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") + scale_fill_gradient2(low = "white", high = "darkblue", mid = "blue", midpoint = 2500, limit = c(0,5000), space = "Lab", name="Dendrogram\nDistance") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()


reorder_cormat <- function(x){
# Use correlation between variables as distance
dd <- as.dist(x)
hc <- hclust(dd)
x <-x[hc$order, hc$order]
}

mphylo <- reorder_cormat(phylo[[4]])
mphylo <- melt(mphylo, na.rm = TRUE)
mphylo <- mphylo[grep("zebrafish", mphylo$Var1),]
mphylo <- mphylo[grep("astyanax", mphylo$Var2),]
ggplot(data = mphylo, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") + scale_fill_gradient2(low = "red", high = "lightskyblue1", mid = "white", midpoint = 500, limit = c(0,1500), space = "Lab", name="Dendrogram\nDistance") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1)) + coord_fixed()


