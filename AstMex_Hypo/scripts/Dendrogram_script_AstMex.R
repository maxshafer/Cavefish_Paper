library(Seurat)
library(dendextend)
library(dplyr)

# Load seurat object (astmex)

load("/Volumes/Maxwell/R_Projects/AstMex_Hypo/AstMex_64k.Robj")

hypo@meta.data$morph_subcluster <- paste(hypo@meta.data$species_morph, hypo@meta.data$SubclusterType)
hypo <- SetAllIdent(hypo, id = "morph_subcluster")
hypo <- BuildClusterTree(hypo)
hypo@cluster.tree[[6]] <- hypo@cluster.tree[[1]]

hypo@meta.data$morph_subtype <- paste(hypo@meta.data$species_morph, hypo@meta.data$Subtype)
hypo <- SetAllIdent(hypo, id = "morph_subtype")
hypo <- BuildClusterTree(hypo)
hypo@cluster.tree[[5]] <- hypo@cluster.tree[[1]]

hypo@meta.data$species_subcluster <- paste(hypo@meta.data$species, hypo@meta.data$SubclusterType)
hypo <- SetAllIdent(hypo, id = "species_subcluster")
hypo <- BuildClusterTree(hypo)
hypo@cluster.tree[[4]] <- hypo@cluster.tree[[1]]

hypo@meta.data$species_subtype <- paste(hypo@meta.data$species, hypo@meta.data$Subtype)
hypo <- SetAllIdent(hypo, id = "species_subtype")
hypo <- BuildClusterTree(hypo)
hypo@cluster.tree[[3]] <- hypo@cluster.tree[[1]]

hypo <- SetAllIdent(hypo, id = "SubclusterType")
hypo <- BuildClusterTree(hypo)
hypo@cluster.tree[[2]] <- hypo@cluster.tree[[1]]

hypo <- SetAllIdent(hypo, id = "Subtype")
hypo <- BuildClusterTree(hypo)


names(hypo@cluster.tree) <- c("Subtype", "SubclusterType", "species_subtype", "species_subcluster", "morph_subtype", "morph_subcluster")

save(hypo, file = "/Volumes/Maxwell/R_Projects/AstMex_Hypo/AstMex_64k.Robj")
