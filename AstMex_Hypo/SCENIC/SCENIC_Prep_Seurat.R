library(reshape2)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(scales)
library(tictoc)

hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")

Idents(hypo.ast) <- "species"
hypo.ast.cave <- subset(hypo.ast, idents = "astyanax_cave")
hypo.ast.surface <- subset(hypo.ast, idents = "astyanax_surface")

# Make a table to ID the min cells from each integrated subcluster (across species)
table <- as.data.frame.matrix(table(hypo.ast@meta.data$Subcluster, hypo.ast@meta.data$species))

table$min <- apply(table, 1, function(x) min(x))

# Generate a list of cell IDs for subsetting

Idents(hypo.ast.cave) <- "Subcluster"
Idents(hypo.ast.surface) <- "Subcluster"

# Find morph specific clusters
prop.table <- table(hypo.ast@meta.data$Subcluster, hypo.ast@meta.data$species)
#prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(hypo.ast@meta.data$Subcluster)
index <- subclusters[!(subclusters %in% c(surface.names, cave.names))]

# These are the same length (which is what we want!)
surface.names <- unlist(lapply(index, function(x) WhichCells(hypo.ast.surface, idents = x, downsample = table[x,3])))
cave.names <- unlist(lapply(index, function(x) WhichCells(hypo.ast.cave, idents = x, downsample = table[x,3])))

# Subset and save for SCENIC analysis
hypo.ast.cave.scenic <- subset(hypo.ast, cells = surface.names)
hypo.ast.surface.scenic <- subset(hypo.ast, cells = cave.names)

saveRDS(hypo.ast.cave.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_cave_20k_vSCENIC.rds")
saveRDS(hypo.ast.surface.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_surface_20k_vSCENIC.rds")

#####################################################
## Prepare the same but for Pachon Tinaja, and Molino

Idents(hypo.ast) <- "morph"
hypo.ast.pachon <- subset(hypo.ast, idents = "Pachon_cave")
hypo.ast.tinaja <- subset(hypo.ast, idents = "Tinaja_cave")
hypo.ast.molino <- subset(hypo.ast, idents = "Molino_cave")

# Make a table to ID the min cells from each integrated subcluster (across species)
table <- as.data.frame.matrix(table(hypo.ast@meta.data$Subcluster, hypo.ast@meta.data$morph))

table$min <- apply(table, 1, function(x) min(x))

# Generate a list of cell IDs for subsetting

Idents(hypo.ast.pachon) <- "Subcluster"
Idents(hypo.ast.tinaja) <- "Subcluster"
Idents(hypo.ast.molino) <- "Subcluster"

# Find morph specific clusters
prop.table <- table(hypo.ast@meta.data$Subcluster, hypo.ast@meta.data$morph)

pachon.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[3]/sum(x) > .9))]
tinaja.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[4]/sum(x) > .9))]
molino.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]

names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2] < 3 | x[3] < 3 | x[4] < 3))]

# Make index
subclusters <- levels(hypo.ast@meta.data$Subcluster)
index <- subclusters[!(subclusters %in% c(names))]

# These are the same length (which is what we want!)
pachon.names <- unlist(lapply(index, function(x) WhichCells(hypo.ast.pachon, idents = x, downsample = table[x,5])))
tinaja.names <- unlist(lapply(index, function(x) WhichCells(hypo.ast.tinaja, idents = x, downsample = table[x,5])))
molino.names <- unlist(lapply(index, function(x) WhichCells(hypo.ast.molino, idents = x, downsample = table[x,5])))

# Subset and save for SCENIC analysis
hypo.ast.pachon.scenic <- subset(hypo.ast, cells = pachon.names)
hypo.ast.tinaja.scenic <- subset(hypo.ast, cells = tinaja.names)
hypo.ast.molino.scenic <- subset(hypo.ast, cells = molino.names)

saveRDS(hypo.ast.pachon.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_pachon_5k_vSCENIC.rds")
saveRDS(hypo.ast.tinaja.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_tinaja_5k_vSCENIC.rds")
saveRDS(hypo.ast.molino.scenic, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_molino_5k_vSCENIC.rds")

