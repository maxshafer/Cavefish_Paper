library(dplyr)
library(tidyr)
library(Seurat)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

## Load marker gene lists

gene.lists <- readRDS(file = "drift_gene_lists_2.rds")

## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4,7,12)] <- lapply(gene.lists[c(1,4,7,12)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]][,2] > 0 & x[[y]][,6] > 0,])) # These elements are from FindConservedMarkers and have different columns
gene.lists.pos[c(2,3,5,6,8,9,10,11,13,14)] <- lapply(gene.lists[c(2,3,5,6,8,9,10,11,13,14)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)

for(i in 1:length(gene.lists.pos)){
  names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

## Filter gene lists by trinarized genes

a = 1.5
b = 2
f = 0.1
trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

## If I change the names of trinarized.exp to match gene.lists.pos (easier to keep rest of script)
## then I can lapply along the trinarized list names, and do another imbedded lapply like below (seeking along cluster IDs)
## Currently doesn't work until I rerun find markers (or rename clusters) - should update names, so I don't have to change this

names <- names(trinarized.exp)

## This now works

gene.lists.pos.2 <- lapply(names, function(y) lapply(names(trinarized.exp[[y]]), function(x) gene.lists.pos[[y]][[x]][rownames(gene.lists.pos[[y]][[x]]) %in% names(trinarized.exp[[y]][[x]]),]))
names(gene.lists.pos.2) <- names(trinarized.exp)

for(i in names){
  names(gene.lists.pos.2[[i]]) <- names(trinarized.exp[[i]])
}

gene.lists.pos.2$cluster.conserved <- lapply( names(gene.lists.pos$cluster.conserved), function(x) gene.lists.pos$cluster.conserved[[x]][rownames(gene.lists.pos$cluster.conserved[[x]]) %in% rownames(gene.lists.pos.2$cluster.zebrafish[[x]]) & rownames(gene.lists.pos$cluster.conserved[[x]]) %in% rownames(gene.lists.pos.2$cluster.astyanax[[x]]),] )
names(gene.lists.pos.2$cluster.conserved) <- names(gene.lists.pos.2$cluster.astyanax)

gene.lists.pos.2$cluster.conserved.ast <- lapply( names(gene.lists.pos$cluster.conserved.ast), function(x) gene.lists.pos$cluster.conserved.ast[[x]][rownames(gene.lists.pos$cluster.conserved.ast[[x]]) %in% rownames(gene.lists.pos.2$cluster.surface[[x]]) & rownames(gene.lists.pos$cluster.conserved.ast[[x]]) %in% rownames(gene.lists.pos.2$cluster.cave[[x]]),] )
names(gene.lists.pos.2$cluster.conserved.ast) <- names(gene.lists.pos.2$cluster.cave)

gene.lists.pos.2$subcluster.conserved <- lapply( names(gene.lists.pos$subcluster.conserved), function(x) gene.lists.pos$subcluster.conserved[[x]][rownames(gene.lists.pos$subcluster.conserved[[x]]) %in% rownames(gene.lists.pos.2$subcluster.zebrafish[[x]]) & rownames(gene.lists.pos$subcluster.conserved[[x]]) %in% rownames(gene.lists.pos.2$subcluster.astyanax[[x]]),] )
names(gene.lists.pos.2$subcluster.conserved) <- names(gene.lists.pos.2$subcluster.astyanax)

gene.lists.pos.2$subcluster.conserved.ast <- lapply( names(gene.lists.pos$subcluster.conserved.ast), function(x) gene.lists.pos$subcluster.conserved.ast[[x]][rownames(gene.lists.pos$subcluster.conserved.ast[[x]]) %in% rownames(gene.lists.pos.2$subcluster.surface[[x]]) & rownames(gene.lists.pos$subcluster.conserved.ast[[x]]) %in% rownames(gene.lists.pos.2$subcluster.cave[[x]]),] )
names(gene.lists.pos.2$subcluster.conserved.ast) <- names(gene.lists.pos.2$subcluster.cave)

gene.lists.pos.2 <- gene.lists.pos.2[names(gene.lists.pos)]

## Save new rds file
saveRDS(gene.lists.pos.2, file = paste("drift_gene_lists_pos_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))
