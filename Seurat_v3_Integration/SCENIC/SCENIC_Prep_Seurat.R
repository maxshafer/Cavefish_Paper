library(reshape2)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(scales)
library(tictoc)

hypo.integrated <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Hypo_integrated_127k_1500VFs_100Dims_vR.rds")
hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vR.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")

# Add integrated labels to hypo.zeb and hypo.ast
# This might introduce NAs (for cells which were removed from hypo.integrated) - these cells shouldn't be included anyway

hypo.zeb@meta.data$integrated_Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(rownames(hypo.zeb@meta.data), rownames(hypo.integrated@meta.data))]
hypo.zeb@meta.data$integrated_Subcluster <- hypo.integrated@meta.data$integrated_Subcluster[match(rownames(hypo.zeb@meta.data), rownames(hypo.integrated@meta.data))]

hypo.ast@meta.data$integrated_Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(rownames(hypo.ast@meta.data), rownames(hypo.integrated@meta.data))]
hypo.ast@meta.data$integrated_Subcluster <- hypo.integrated@meta.data$integrated_Subcluster[match(rownames(hypo.ast@meta.data), rownames(hypo.integrated@meta.data))]


# Make a table to ID the min cells from each integrated subcluster (across species)
table <- as.data.frame.matrix(table(hypo.integrated@meta.data$integrated_Subcluster, hypo.integrated@meta.data$species.2))

table$min <- apply(table, 1, function(x) min(x))

# Generate a list of cell IDs for subsetting

Idents(hypo.ast) <- "integrated_Subcluster"
Idents(hypo.zeb) <- "integrated_Subcluster"

# Find species specific clusters
prop.table <- table(hypo.integrated@meta.data$integrated_Subcluster, hypo.integrated@meta.data$species.2)
#prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
zeb.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
ast.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo.integrated))
index <- subclusters[!(subclusters %in% c(zeb.names, ast.names))]

# These are the same length (which is what we want!)
ast.names <- unlist(lapply(index, function(x) WhichCells(hypo.ast, idents = x, downsample = table[x,3])))
zeb.names <- unlist(lapply(index, function(x) WhichCells(hypo.zeb, idents = x, downsample = table[x,3])))

# Subset and save for SCENIC analysis
hypo.zeb.scenic <- subset(hypo.zeb, cells = zeb.names)
hypo.ast.scenic <- subset(hypo.ast, cells = ast.names)

saveRDS(hypo.zeb.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vSCENIC.rds")
saveRDS(hypo.ast.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vSCENIC.rds")

