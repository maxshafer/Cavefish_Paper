library(Seurat)

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


conserved.markers <- lapply(levels(Idents(hypo.integrated.zeb))[!(levels(Idents(hypo.integrated.zeb)) %in% "Cilliated")], function(x) FindConservedMarkers(hypo.integrated, ident.1 = x, grouping.var = "species.2", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers) <- levels(Idents(hypo.integrated.zeb))[!(levels(Idents(hypo.integrated.zeb)) %in% "Cilliated")]

conserved.markers.ast <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 1000))
names(conserved.markers.ast) <- levels(Idents(hypo.integrated.ast))

# Find Markers for other subsets

zebrafish.markers <- lapply(levels(Idents(hypo.integrated.zeb))[!(levels(Idents(hypo.integrated.zeb)) %in% "Cilliated")], function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(zebrafish.markers) <- levels(Idents(hypo.integrated.zeb))[!(levels(Idents(hypo.integrated.zeb)) %in% "Cilliated")]

astyanax.markers <- lapply(levels(Idents(hypo.integrated.ast)), function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(astyanax.markers) <- levels(Idents(hypo.integrated.ast))
astyanax.markers <- astyanax.markers[c(names(zebrafish.markers))]

surface.markers <- lapply(levels(Idents(hypo.integrated.surface)), function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(surface.markers) <- levels(Idents(hypo.integrated.surface))

cave.markers <- lapply(levels(Idents(hypo.integrated.cave)), function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = T, max.cells.per.ident = 1000))
names(cave.markers) <- levels(Idents(hypo.integrated.cave))

gene.lists <- list(conserved.markers, zebrafish.markers, astyanax.markers, conserved.markers.ast, surface.markers, cave.markers)

saveRDS(gene.lists, file = "drift_gene_lists.rds")

############## SubclusterType Markers ##############

## Find Conserved Sub Markers
## First identify those subclusters that are species specific

## NEED TO REMAKE INDEX FILE FOR ZEB AST COMPARISON

Idents(hypo.integrated) <- "integrated_SubclusterType"
Idents(hypo.integrated.zeb) <- "integrated_SubclusterType"
Idents(hypo.integrated.ast) <- "integrated_SubclusterType"

# Find species specific clusters
prop.table <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species.2)
#prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
zeb.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
ast.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo.integrated))
index <- subclusters[!(subclusters %in% c(zeb.names, ast.names))]

# Find Conserved Sub markers across species.2

conserved.markers.sub <- lapply(index, function(x) FindConservedMarkers(hypo.integrated, ident.1 = x, grouping.var = "species.2", verbose = T, max.cells.per.ident = 500))
names(conserved.markers.sub) <- index

gene.lists[[7]] <- conserved.markers.sub

# Find Sub Markers for each species.2

zebrafish.markers.sub <- lapply(index, function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(zebrafish.markers.sub) <- index

zeb.name.markers <- lapply(zeb.names, function(x) FindMarkers(hypo.integrated.zeb, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(zeb.name.markers) <- zeb.names

gene.lists[[8]] <- zebrafish.markers.sub

astyanax.markers.sub <- lapply(index, function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(astyanax.markers.sub) <- index

ast.name.markers <- lapply(ast.names, function(x) FindMarkers(hypo.integrated.ast, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(ast.name.markers) <- ast.names

gene.lists[[9]] <- astyanax.markers.sub

gene.lists[[10]] <- zeb.name.markers
gene.lists[[11]] <- ast.name.markers

saveRDS(gene.lists, file = "drift_gene_lists.rds")

### Find Sub Markers for each species (cave/surface)

gene.lists <- readRDS("drift_gene_lists.rds")

Idents(hypo.integrated.ast) <- "integrated_SubclusterType"
Idents(hypo.integrated.surface) <- "integrated_SubclusterType"
Idents(hypo.integrated.cave) <- "integrated_SubclusterType"

# Find morph specific clusters
prop.table <- table(hypo.integrated.ast@meta.data$integrated_SubclusterType, hypo.integrated.ast@meta.data$species)
# prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo.integrated.ast))
index2 <- subclusters[!(subclusters %in% c(zeb.names, surface.names, cave.names))]

conserved.markers.ast.sub <- lapply(index2, function(x) FindConservedMarkers(hypo.integrated.ast, ident.1 = x, grouping.var = "species", verbose = T, max.cells.per.ident = 500))
names(conserved.markers.ast.sub) <- index2

gene.lists[[12]] <- conserved.markers.ast.sub

surface.markers.sub <- lapply(index2, function(x) FindMarkers(hypo.integrated.surface, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(surface.markers.sub) <- index2

gene.lists[[13]] <- surface.markers.sub

cave.markers.sub <- lapply(index2, function(x) FindMarkers(hypo.integrated.cave, ident.1 = x, verbose = T, max.cells.per.ident = 500))
names(cave.markers.sub) <- index2

gene.lists[[14]] <- cave.markers.sub

names(gene.lists) <- c("conserved.markers", "zebrafish.markers", "astyanax.markers", "conserved.markers.ast", "surface.markers", "cave.markers", "conserved.markers.sub", "zebrafish.markers.sub", "astyanax.markers.sub", "zeb.name.markers", "ast.name.markers", "conserved.markers.ast.sub", "surface.markers.sub", "cave.markers.sub")

saveRDS(gene.lists, file = "drift_gene_lists.rds")


## Remove weird subclusters (from species-specific analysis)
## These were likely created by the integration process (hbaa2 and immune/blood genes in neuronal clusters)
## Save as V2 (change old V2 as V1, and old V1 as V0)

Idents(hypo.integrated) <- "integrated_SubclusterType"

hypo.integrated <- subset(hypo.integrated, cells = WhichCells(hypo.integrated, idents = c("Glut_3_5", "GABA_1_12", "Glut_3_7", "Glut_5_5"), invert = T))
hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_3_6"] <- "Glut_3_5"





hypo.integrated@meta.data$integrated_SubclusterType <- factor(hypo.integrated@meta.data$integrated_SubclusterType, levels(hypo.integrated@meta.data$integrated_SubclusterType)[!(levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% c("Glut_3_6", "GABA_1_12", "Glut_3_7", "Glut_5_5"))])

saveRDS(hypo.integrated, file = "Hypo_integrated_127k_1500VFs_100Dims_v2.rds")

## Combine Glut_2 clusters into 3 clusters
## Glut_2_0 + Glut_2_2, Glut_2_3 + Glut_2_4, and Glut_2_1, and Glut_2_5
## Find new marker genes for Glut_2_0, Glut_2_2, Glut_2_3 and replace in gene.lists

Idents(hypo.integrated) <- "integrated_SubclusterType"

hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_2_2"] <- "Glut_2_0"
hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_2_3"] <- "Glut_2_2"
hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_2_4"] <- "Glut_2_2"

hypo.integrated@meta.data$integrated_SubclusterType[hypo.integrated@meta.data$integrated_SubclusterType == "Glut_2_5"] <- "Glut_2_3"

hypo.integrated@meta.data$integrated_SubclusterType <- factor(hypo.integrated@meta.data$integrated_SubclusterType, levels(hypo.integrated@meta.data$integrated_SubclusterType)[!(levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% c("Glut_2_4", "Glut_2_5"))])

gene.lists[["conserved.markers.sub"]][["Glut_2_2"]] <- FindConservedMarkers(object = hypo.integrated, ident.1 = "Glut_2_2", grouping.var = "species.2", verbose = T, max.cells.per.ident = 500)
gene.lists[["zebrafish.markers.sub"]][["Glut_2_2"]] <- FindMarkers(hypo.integrated.zeb, ident.1 = "Glut_2_2", verbose = T, max.cells.per.ident = 500)
gene.lists[["astyanax.markers.sub"]][["Glut_2_2"]] <- FindMarkers(hypo.integrated.ast, ident.1 = "Glut_2_2", verbose = T, max.cells.per.ident = 500)
gene.lists[["conserved.markers.sub"]][["Glut_2_3"]] <- FindConservedMarkers(object = hypo.integrated, ident.1 = "Glut_2_3", grouping.var = "species.2", verbose = T, max.cells.per.ident = 500)
gene.lists[["zebrafish.markers.sub"]][["Glut_2_3"]] <- FindMarkers(hypo.integrated.zeb, ident.1 = "Glut_2_3", verbose = T, max.cells.per.ident = 500)
gene.lists[["astyanax.markers.sub"]][["Glut_2_3"]] <- FindMarkers(hypo.integrated.ast, ident.1 = "Glut_2_3", verbose = T, max.cells.per.ident = 500)
gene.lists[["zeb.name.markers"]][["Glut_2_0"]] <- FindMarkers(hypo.integrated.zeb, ident.1 = "Glut_2_0", verbose = T, max.cells.per.ident = 500)

gene.lists[["conserved.markers.sub"]] <- gene.lists[["conserved.markers.sub"]][levels(hypo.integrated@meta.data$integrated_SubclusterType)[levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% names(gene.lists[["conserved.markers.sub"]])]]
gene.lists[["zebrafish.markers.sub"]] <- gene.lists[["zebrafish.markers.sub"]][levels(hypo.integrated@meta.data$integrated_SubclusterType)[levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% names(gene.lists[["zebrafish.markers.sub"]])]]
gene.lists[["astyanax.markers.sub"]] <- gene.lists[["astyanax.markers.sub"]][levels(hypo.integrated@meta.data$integrated_SubclusterType)[levels(hypo.integrated@meta.data$integrated_SubclusterType) %in% names(gene.lists[["astyanax.markers.sub"]])]]
gene.lists[["zeb.name.markers"]] <- gene.lists[["zeb.name.markers"]][names(gene.lists[["zeb.name.markers"]]) %in% zeb.names]
gene.lists[["ast.name.markers"]] <- gene.lists[["ast.name.markers"]][names(gene.lists[["ast.name.markers"]]) %in% ast.names]

saveRDS(gene.lists, file = "drift_gene_lists.rds")

saveRDS(hypo.integrated, file = "Hypo_integrated_127k_1500VFs_100Dims_v3.rds")
