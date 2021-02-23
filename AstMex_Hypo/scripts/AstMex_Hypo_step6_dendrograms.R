library(Seurat)
library(Matrix)

# Load seurat object (astmex)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo <- readRDS("AstMex_63k_vR.rds")
DefaultAssay(hypo) <- "RNA"


# Find morph specific clusters
prop.table <- table(hypo@meta.data$Subcluster, hypo@meta.data$species)
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(hypo@meta.data$Subcluster)
index <- subclusters[!(subclusters %in% c(surface.names, cave.names))]

# Make Dendrograms as a new list

dendrograms <- list()

Idents(hypo) <- "Cluster"
hypo <- BuildClusterTree(hypo)
dendrograms[[1]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "Subcluster"
hypo <- BuildClusterTree(hypo)
dendrograms[[2]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "species_Cluster"
hypo <- BuildClusterTree(hypo)
dendrograms[[3]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "Subcluster"
hypo.1 <- subset(hypo, idents = index)
Idents(hypo.1) <- "species_Subcluster"
hypo.1 <- BuildClusterTree(hypo.1)
dendrograms[[4]] <- Tool(hypo.1, slot = "BuildClusterTree")

Idents(hypo) <- "morph_Cluster"
hypo <- BuildClusterTree(hypo)
dendrograms[[5]] <- Tool(hypo, slot = "BuildClusterTree")

Idents(hypo) <- "morph_Subcluster"
hypo <- BuildClusterTree(hypo)
dendrograms[[6]] <- Tool(hypo, slot = "BuildClusterTree")

names(dendrograms) <- c("Cluster", "Subcluster", "species_Cluster", "species_Subcluster", "morph_Cluster", "morph_Subcluster")

saveRDS(dendrograms, file = "Ast_dendrograms.rds")
