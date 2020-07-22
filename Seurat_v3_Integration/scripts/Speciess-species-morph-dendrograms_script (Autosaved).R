library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)
library(png)
library(grid)
library(ggtree)
library(magick)
library(ggimage)
library(phytools)
library(ggplotify)


# load objects

load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj")
load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_64k.Robj")

genes <- intersect(row.names(hypo.ast@data), row.names(hypo.zeb@data))
# OR
# genes2 <- union(hypo.ast@var.genes, hypo.zeb@var.genes)
# genes2 <- genes2[genes2 %in% row.names(hypo.ast@data)]
# genes2 <- genes2[genes2 %in% row.names(hypo.zeb@data)]

# genes <- genes2

# Make density phylogram

hypo.combined <- MergeSeurat(hypo.zeb, hypo.ast)
hypo.combined@meta.data$species_morph[hypo.combined@meta.data$species == "zebrafish"] <- "zebrafish"


hypo.combined <- SetAllIdent(hypo.combined, id = "species_morph")


## For each of idents, need to subset and BuildClusterTree

type2 <- lapply(unique(hypo.combined@meta.data$Type2), function(x) BuildClusterTree(SubsetData(hypo.combined, cells.use = WhichCells(hypo.combined, subset.name = "Type2", accept.value = x)), genes.use = genes))
type2 <- lapply(type2, function(x) x@cluster.tree[[1]])

names(type2) <- unique(hypo.combined@meta.data$Type2)
class(type2) <- "multiPhylo"
#ggtree(type2) + facet_wrap(~.id, scale = "free") + geom_tiplab()

densityTree(type2, type = "phylogram", use.edge.length = F, nodes = "inner", show.axis = F)



## For each of idents, need to subset and BuildClusterTree

subtype <- lapply(unique(hypo.combined@meta.data$integrated_Subtype), function(x) BuildClusterTree(SubsetData(hypo.combined, cells.use = WhichCells(hypo.combined, subset.name = "integrated_Subtype", accept.value = x)), genes.use = genes))
subtype <- lapply(subtype, function(x) x@cluster.tree[[1]])
names(subtype) <- idents
class(subtype) <- "multiPhylo"
ggtree(subtype) + facet_wrap(~.id, scale = "free") + geom_tiplab()

densityTree(subtype, type = "phylogram", use.edge.length = F, nodes = "inner", show.axis = F)




## For subclusters

SubclusterTypes <- unique(hypo.combined@meta.data$integrated_SubclusterType)
index <- table(hypo.combined@meta.data$integrated_SubclusterType, hypo.combined@meta.data$species_morph)
index <- t(apply(index, 1, function(x) x > 2))
index <- apply(index, 1, function(x) length(x[x]) > 2)

length(SubclusterTypes[SubclusterTypes %in% names(index)[index]])


## use group_by on above to actually count the number of cells, and skip the below?

subcluster <- lapply(SubclusterTypes[SubclusterTypes %in% names(index)[index]], function(x) BuildClusterTree(SubsetData(hypo.combined, cells.use = WhichCells(hypo.combined, subset.name = "integrated_SubclusterType", accept.value = x)), genes.use = genes))


subcluster <- lapply(subcluster, function(x) x@cluster.tree[[1]])
names(subcluster) <- SubclusterTypes[SubclusterTypes %in% names(index)[index]]
class(subcluster) <- "multiPhylo"
ggtree(subcluster) + facet_wrap(~.id, scale = "free") + geom_tiplab()

ggdensitree(subcluster, alpha = 0.3, colour = "black", layout = "rectangular", align.tips = T, branch.length = "none") + geom_tiplab(size = 3) + theme(plot.margin = margin(1,1,1,1, "cm"))


densityTree(subcluster, type = "cladogram", use.edge.length = F, nodes = "inner", show.axis = F, alpha = 0.03)

test <- c(type2, subtype, subcluster)
densityTree(test, type = "cladogram", use.edge.length = F, nodes = "inner", show.axis = F, alpha = 0.05)


