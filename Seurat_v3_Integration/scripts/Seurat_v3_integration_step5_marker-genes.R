library(Seurat)
library(Matrix)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_v1.rds")

DefaultAssay(hypo.integrated) <- "RNA"

Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

Idents(hypo.integrated.ast) <- "species"
hypo.integrated.surface <- subset(hypo.integrated.ast, idents = "astyanax_surface")
hypo.integrated.cave <- subset(hypo.integrated.ast, idents = "astyanax_cave")

############## Subtype Markers ##############

# Find Conserved Markers

Idents(hypo.integrated) <- "integrated_Subtype"
Idents(hypo.integrated.zeb) <- "integrated_Subtype"
Idents(hypo.integrated.ast) <- "integrated_Subtype"
Idents(hypo.integrated.surface) <- "integrated_Subtype"
Idents(hypo.integrated.cave) <- "integrated_Subtype"


conserved.markers <- lapply(levels(Idents(hypo.integrated)), function(x) FindConservedMarkers(hypo.integrated, ident.1 = x, grouping.var = "species.2", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers) <- levels(Idents(hypo.integrated))

conserved.markers.ast <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers.ast) <- levels(Idents(hypo.integrated.ast))

# Find Markers for other subsets

zebrafish.markers <- lapply(levels(Idents(hypo.integrated.zeb)), function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(zebrafish.markers) <- levels(Idents(hypo.integrated.zeb))

astyanax.markers <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(astyanax.markers) <- levels(Idents(hypo.integrated.ast))
astyanax.markers <- astyanax.markers[c(names(zebrafish.markers))]

surface.markers <- lapply(levels(Idents(hypo.integrated.surface)), function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(surface.markers) <- levels(Idents(hypo.integrated.surface))

cave.markers <- lapply(levels(Idents(hypo.integrated.cave)), function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(cave.markers) <- levels(Idents(hypo.integrated.cave))

gene.lists <- list(conserved.markers, zebrafish.markers, astyanax.markers, surface.markers, cave.markers)

saveRDS(gene.lists, file = "drift_gene_lists.rds")

############## SubclusterType Markers ##############

## Find Conserved Sub Markers
## First identify those subclusters that are species specific

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

Idents(hypo.integrated) <- "integrated_SubclusterType"
Idents(hypo.integrated.zeb) <- "integrated_SubclusterType"
Idents(hypo.integrated.ast) <- "integrated_SubclusterType"

# Find morph specific clusters
prop.table <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species.2)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$zeb_specific <- ifelse(prop.table$zebrafish > .9, "yes", "no")
prop.table$ast_specific <- ifelse(prop.table$astyanax > .9, "yes", "no")

zeb.names <- row.names(prop.table[prop.table$zeb_specific == "yes",])
ast.names <- row.names(prop.table[prop.table$ast_specific == "yes",])

# Make index
subclusters <- levels(Idents(hypo.integrated))
index <- subclusters[!(subclusters %in% c(zeb.names, ast.names))]

# Find Conserved Sub markers across species.2

conserved.markers.sub <- lapply(index, function(x) FindConservedMarkers(hypo.integrated, ident.1 = x, grouping.var = "species.2", verbose = T, max.cells.per.ident = 500))
names(conserved.markers.sub) <- index

gene.lists[[6]] <- conserved.markers.sub

# Find Sub Markers for each species.2

zebrafish.markers.sub <- lapply(index, function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(zebrafish.markers.sub) <- index

zeb.name.markers <- lapply(zeb.names, function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(zeb.name.markers) <- zeb.names

gene.lists[[7]] <- zebrafish.markers.sub

astyanax.markers.sub <- lapply(index, function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(astyanax.markers.sub) <- index

ast.name.markers <- lapply(ast.names, function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(ast.name.markers) <- ast.names

gene.lists[[8]] <- astyanax.markers.sub

saveRDS(gene.lists, file = "drift_gene_lists.rds")

### Find Sub Markers for each species (cave/surface)

gene.lists <- readRDS("drift_gene_lists.rds")

Idents(hypo.integrated.ast) <- "integrated_SubclusterType"
Idents(hypo.integrated.surface) <- "integrated_SubclusterType"
Idents(hypo.integrated.cave) <- "integrated_SubclusterType"

# Find morph specific clusters
prop.table <- table(hypo.integrated.ast@meta.data$integrated_SubclusterType, hypo.integrated.ast@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$surface_specific <- ifelse(prop.table$astyanax_surface > .9, "yes", "no")
prop.table$cave_specific <- ifelse(prop.table$astyanax_cave > .9, "yes", "no")

# prop.table <- prop.table[index,]

surface.names <- row.names(prop.table[prop.table$surface_specific == "yes",])
cave.names <- row.names(prop.table[prop.table$cave_specific == "yes",])

# Make index
subclusters <- levels(Idents(hypo.integrated.ast))
index2 <- subclusters[!(subclusters %in% c(zeb.names, surface.names, cave.names, "Oligodendrocytes_5"))]

conserved.markers.ast.sub <- lapply(index2, function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 500))
names(conserved.markers.ast.sub) <- index2

gene.lists[[9]] <- conserved.markers.ast.sub

surface.markers.sub <- lapply(index2, function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(surface.markers.sub) <- index2

gene.lists[[10]] <- surface.markers.sub

cave.markers.sub <- lapply(index2, function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(cave.markers.sub) <- index2

gene.lists[[11]] <- cave.markers.sub

saveRDS(gene.lists, file = "drift_gene_lists.rds")

# gene.lists <- list(conserved.markers, zebrafish.markers, astyanax.markers, conserved.markers.sub, zebrafish.markers.sub, astyanax.markers.sub, conserved.markers.ast, surface.markers, cave.markers, conserved.markers.ast.sub, surface.markers.sub, cave.markers.sub)


DotPlot(hypo.integrated, features = head(row.names(zeb.name.markers[["Glut_1_3"]]), 10), group.by = "integrated_SubclusterType") + RotatedAxis() + coord_flip() + theme(axis.text = element_text(size = 8))

DotPlot(hypo.integrated, features = head(row.names(ast.name.markers[["Progenitors_11"]]), 10), group.by = "integrated_SubclusterType") + RotatedAxis() + coord_flip() + theme(axis.text = element_text(size = 8))

DotPlot(hypo.integrated, features = c("ENSAMXG00000025407", "ENSAMXG00000007964", "ENSAMXG00000015001"), group.by = "integrated_SubclusterType") + RotatedAxis() + coord_flip() + theme(axis.text = element_text(size = 8))

DotPlot(hypo.integrated, features = head(row.names(gene.lists[[1]][["Glut_0"]]), 10), group.by = "integrated_SubclusterType") + RotatedAxis() + coord_flip() + theme(axis.text = element_text(size = 8))

