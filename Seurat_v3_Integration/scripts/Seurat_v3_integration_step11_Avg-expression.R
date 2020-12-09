library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

# load objects

hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")
hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v3.rds")


## Add integrated ids to objects

hypo.zeb@meta.data$integrated_Subtype <- hypo.integrated@meta.data$integrated_Subtype[match(row.names(hypo.zeb@meta.data), row.names(hypo.integrated@meta.data))]

hypo.zeb@meta.data$integrated_SubclusterType <- hypo.integrated@meta.data$integrated_SubclusterType[match(row.names(hypo.zeb@meta.data), row.names(hypo.integrated@meta.data))]

hypo.ast@meta.data$integrated_Subtype <- hypo.integrated@meta.data$integrated_Subtype[match(row.names(hypo.ast@meta.data), row.names(hypo.integrated@meta.data))]

hypo.ast@meta.data$integrated_SubclusterType <- hypo.integrated@meta.data$integrated_SubclusterType[match(row.names(hypo.ast@meta.data), row.names(hypo.integrated@meta.data))]

# rm(hypo.integrated)

# saveRDS(hypo.zeb, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
# saveDS(hypo.ast, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")

## Make cave and surface objects, then the list

Idents(hypo.ast) <- "species"

hypo.cave <- subset(hypo.ast, idents = "astyanax_cave")
hypo.surface <- subset(hypo.ast, idents = "astyanax_surface")

hypo.list <- list(hypo.zeb, hypo.ast, hypo.surface, hypo.cave)

# ## To do for each of the 16 zeb samples, add to hypo list, and just run for Subtypes
# 
# Idents(hypo.zeb) <- "orig.ident"
# hypo.samples <- lapply(levels(hypo.zeb@meta.data$orig.ident), function(x) subset(hypo.zeb, idents = x))
# 
# ## To do for each of the ast samples, add to hypo list, and just run for Subtypes/Types
# 
# Idents(hypo.ast) <- "orig.ident"
# hypo.samples2 <- lapply(levels(hypo.ast@meta.data$orig.ident), function(x) subset(hypo.ast, idents = x))
# 
# hypo.list <- c(hypo.list, hypo.samples, hypo.samples2)

## For each object (zeb, ast, surface, cave), and for each category (Type2, Subtype, SubclusterType), calculate averages and make a matrix (list of lists of matrices)

categories <- c("Subtype", "SubclusterType", "integrated_Subtype", "integrated_SubclusterType")

norm.cluster <- list()
for(i in 1:length(categories)) {
	for (j in 1:length(hypo.list)){
	  Idents(hypo.list[[j]]) <- categories[[i]]
	}
	idents.list <- lapply(hypo.list, function(x) levels(Idents(x))) # this is a list of the Type2 identities for each object, not the categories
	idents.list <- lapply(seq_along(idents.list), function(x) idents.list[[x]][as.vector(table(Idents(hypo.list[[x]]))) > 1])
	
	norm.cluster[[i]] <- lapply(seq_along(hypo.list), function(x) lapply(seq_along(idents.list[[x]]), function(y) data.frame(mean.exp = apply(expm1(GetAssayData(hypo.list[[x]])[,WhichCells(hypo.list[[x]], idents = idents.list[[x]][[y]])]), 1, function(z) mean(z)))))

	names(norm.cluster[[i]]) <- c("hypo.zeb", "hypo.ast", "hypo.surface", "hypo.cave")

	for(j in 1:length(hypo.list)) {
		names(norm.cluster[[i]][[j]]) <- idents.list[[j]]
	}
}

names(norm.cluster) <- categories

for(i in 1:length(categories)) {
  for (j in 1:length(hypo.list)){
    Idents(hypo.list[[j]]) <- categories[[i]]
  }
  idents.list <- lapply(hypo.list, function(x) levels(Idents(x))) # this is a list of the Type2 identities for each object, not the categories
  idents.list <- lapply(seq_along(idents.list), function(x) idents.list[[x]][as.vector(table(Idents(hypo.list[[x]]))) > 1])

  for(j in 1:length(hypo.list)) {
    colnames(norm.cluster[[i]][[j]]) <- idents.list[[j]]
  }
}

saveRDS(norm.cluster, file = "Normed_expression_data.rds")

# ## After quality control, several cell types were removed, remove them from this list
# ## These likely represent artifacts from integration: They are likely erthrocytes which were batch correctedwrongly
# normed.expression <- readRDS("Normed_expression_data.rds")
# 
# remove.types <- c("GABA_1_12", "Glut_3_5", "Glut_3_7", "Glut_5_5")
# 
# normed.expression[[4]] <- lapply(normed.expression[[4]], function(x) x[!(names(x) %in% remove.types)])
# 
# saveRDS(normed.expression, file = "Normed_expression_data.rds")
