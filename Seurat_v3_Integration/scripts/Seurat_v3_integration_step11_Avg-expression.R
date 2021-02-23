library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

# load objects

hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vR.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")
hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")


## Add integrated ids to objects

hypo.zeb@meta.data$integrated_Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(row.names(hypo.zeb@meta.data), row.names(hypo.integrated@meta.data))]
hypo.zeb@meta.data$integrated_Subcluster <- hypo.integrated@meta.data$integrated_Subcluster[match(row.names(hypo.zeb@meta.data), row.names(hypo.integrated@meta.data))]
hypo.ast@meta.data$integrated_Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(row.names(hypo.ast@meta.data), row.names(hypo.integrated@meta.data))]
hypo.ast@meta.data$integrated_Subcluster <- hypo.integrated@meta.data$integrated_Subcluster[match(row.names(hypo.ast@meta.data), row.names(hypo.integrated@meta.data))]

# rm(hypo.integrated)

saveRDS(hypo.zeb, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vR.rds")
saveDS(hypo.ast, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")

## Make cave and surface objects, then the list

Idents(hypo.ast) <- "species"

hypo.cave <- subset(hypo.ast, idents = "astyanax_cave")
hypo.surface <- subset(hypo.ast, idents = "astyanax_surface")

hypo.list <- list(hypo.zeb, hypo.ast, hypo.surface, hypo.cave)

# ## To do for each of the 16 zeb samples, add to hypo list, and just run for Clusters
# 
# Idents(hypo.zeb) <- "orig.ident"
# hypo.samples <- lapply(levels(hypo.zeb@meta.data$orig.ident), function(x) subset(hypo.zeb, idents = x))
# 
# ## To do for each of the ast samples, add to hypo list, and just run for Clusters/Types
# 
# Idents(hypo.ast) <- "orig.ident"
# hypo.samples2 <- lapply(levels(hypo.ast@meta.data$orig.ident), function(x) subset(hypo.ast, idents = x))
# 
# hypo.list <- c(hypo.list, hypo.samples, hypo.samples2)

## For each object (zeb, ast, surface, cave), and for each category (Type2, Cluster, Subcluster), calculate averages and make a matrix (list of lists of matrices)

categories <- c("Cluster", "Subcluster", "integrated_Cluster", "integrated_Subcluster")

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

# Need to cbind the lists, then rename the columns? Is that what the below does?


for(i in 1:length(categories)) {
  for (j in 1:length(hypo.list)){
    Idents(hypo.list[[j]]) <- categories[[i]]
  }
  idents.list <- lapply(hypo.list, function(x) levels(Idents(x)))
  idents.list <- lapply(seq_along(idents.list), function(x) idents.list[[x]][as.vector(table(Idents(hypo.list[[x]]))) > 1])

  for(j in 1:length(hypo.list)) {
    norm.cluster[[i]][[j]] <- Reduce(cbind, norm.cluster[[i]][[j]])
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
