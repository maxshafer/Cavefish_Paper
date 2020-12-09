library(dplyr)
library(tidyr)
library(Seurat)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_v3.rds")
Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

## Load marker gene lists

gene.lists <- readRDS(file = "drift_gene_lists.rds")

## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4,7,12)] <- lapply(gene.lists[c(1,4,7,12)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]][,2] > 0 & x[[y]][,6] > 0,])) # These elements are from FindConservedMarkers and have different columns
gene.lists.pos[c(2,3,5,6,8,9,10,11,13,14)] <- lapply(gene.lists[c(2,3,5,6,8,9,10,11,13,14)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists)

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

saveRDS(gene.lists.pos, file = "drift_gene_lists_pos.rds") 

gene.lists.pos <- readRDS("drift_gene_lists_pos.rds")

## Calculate the drift index!
## Takes 3 lists of marker genes, conserved, speices.1 and species.2
## Takes the intersection of all 3 Subtype names, and subset/order lists by that
## Then run Drift Index caculation for each Subtype
## Maybe add a part for using a cutoff (p-val, logfc, percent)
## Add a part for putting number of genes

calcDriftIndex <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	DI <- list()
	if (is.null(subset)) {
		DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (nrow(conserved[[x]]) / nrow(species.1[[x]]))) * (1 - (nrow(conserved[[x]]) / nrow(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
	}
	}
	names(DI) <- names
	return(DI)
}

## For Subtypes
drift.index.Subtypes <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers, species.1 = gene.lists.pos$zebrafish.markers, species.2 = gene.lists.pos$astyanax.markers)
DI <- data.frame(Subtype = names(drift.index.Subtypes), values = unlist(drift.index.Subtypes))

drift.index.Subtypes.ast <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers.ast, species.1 = gene.lists.pos$surface.markers, species.2 = gene.lists.pos$cave.markers)
DI.ast <- data.frame(Subtype = c(names(unlist(drift.index.Subtypes.ast))), values = c(unlist(drift.index.Subtypes.ast)))

## For Subclusters
drift.index.SubclusterType <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub)
DI.sub <- data.frame(SubclusterType = c(names(unlist(drift.index.SubclusterType))), values = c(unlist(drift.index.SubclusterType)))

drift.index.SubclusterTypes.ast <- calcDriftIndex(conserved = gene.lists.pos$conserved.markers.ast.sub, species.1 = gene.lists.pos$surface.markers.sub, species.2 = gene.lists.pos$cave.markers.sub)
DI.ast.sub <- data.frame(SubclusterType = c(names(unlist(drift.index.SubclusterTypes.ast))), values = c(unlist(drift.index.SubclusterTypes.ast)))

# OK, add Subtype as a column to subcluster df

index <- unique(data.frame(Subtype = hypo.integrated@meta.data$integrated_Subtype, SubclusterType = hypo.integrated@meta.data$integrated_SubclusterType))
DI.sub$Subtype <- index$Subtype[match(DI.sub$SubclusterType, index$SubclusterType)]
DI.sub$Subtype <- factor(DI.sub$Subtype, levels = levels(hypo.integrated@meta.data$integrated_Subtype))

index.ast <- unique(tibble(Subtype = hypo.integrated.ast@meta.data$integrated_Subtype, SubclusterType = hypo.integrated.ast@meta.data$integrated_SubclusterType))
DI.ast.sub$Subtype <- index.ast$Subtype[match(DI.ast.sub$SubclusterType, index.ast$SubclusterType)]
DI.ast.sub$Subtype <- factor(DI.ast.sub$Subtype, levels = levels(hypo.integrated@meta.data$integrated_Subtype))

# Caculate drift for only TFs, NPS/NTS etc
# Get Go lists and use only those that are marker genes

go_lists <- list.files("../Seurat_v3_Integration/SCENIC")[grep("GO", list.files("../Seurat_v3_Integration/SCENIC"))]
go_lists <- lapply(paste("../Seurat_v3_Integration/SCENIC/", go_lists, sep = ""), function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files("../Seurat_v3_Integration/SCENIC")[grep("GO", list.files("../Seurat_v3_Integration/SCENIC"))]
names(go_lists)
go_lists[[13]] <- Reduce(union, list(go_lists[[5]], go_lists[[6]], go_lists[[9]]))

## Plot Subtypes

DI2 <- tibble(cell_type = c(names(gene.lists.pos$conserved.markers.sub)), all = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, invert = F))), TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = go_lists[[3]], invert = F))), NP_NTS = c(unlist(calcDriftIndex(conserved = gene.lists.pos$conserved.markers.sub, species.1 = gene.lists.pos$zebrafish.markers.sub, species.2 = gene.lists.pos$astyanax.markers.sub, subset = go_lists[[13]], invert = F))))
DI2 <- DI2[!is.na(DI2$NP_NTS),]
DI.sub.GO <- reshape2::melt(DI2)
DI.sub.GO$variable <- factor(DI.sub.GO$variable, levels = c("all", "NP_NTS", "TFs"))
DI.sub.GO$value[is.infinite(DI.sub.GO$value)] <- 0

## Save

DI.list <- list(DI, DI.sub, DI.ast, DI.ast.sub, DI.sub.GO)
names(DI.list) <- c("DI", "DI.sub", "DI.ast", "DI.ast.sub", "DI.sub.GO")

saveRDS(DI.list, file = "DI_results.rds")


## Calculate Drift between all cell types both within and between species

## OK, conserved marker genes are found by comparing the p vals between lists
## can use metap::minimump to exclude genes which appear in both gene lists, but are not significant
## compares the two p_val (not adj_p_val), and outputs a minimump_p_val
## OR just the intersect, or calculate the intersect, and then also the minimump
## so need to make data.frames for all cell types comparisons? and then compute the minimump

## Populate a matrix with the DI number, using the row and col names

calcDriftIndex2 <- function(cell_type1 = cell_type1, cell_type2 = cell_type2) {
  colnames(cell_type1) <- paste("cell_type1", colnames(cell_type1), sep = "_")
  colnames(cell_type2) <- paste("cell_type2", colnames(cell_type2), sep = "_")
  index <- intersect(row.names(cell_type1), row.names(cell_type2))
  df <- cbind(cell_type1[index,], cell_type2[index,])
  df$minimup_p_val <- apply(df[,grep("p_val$", colnames(df))], 1, function(x) metap::minimump(x)$p)
  df <- df[df$minimup_p_val < 0.05,]
  
  DI <- 1 - sqrt( abs( (1 - (nrow(df) / nrow(cell_type1))) * (1 - (nrow(df) / nrow(cell_type2))) ) )
  return(DI)
}

matrix.zeb <- list()
for (i in 1:length(gene.lists.pos$zebrafish.markers)) {
  matrix.zeb[[i]] <- lapply(gene.lists.pos$zebrafish.markers, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$zebrafish.markers[[i]], cell_type2 = x))
}
DI.list[[6]] <- matrix.zeb

matrix.ast <- list()
for (i in 1:length(gene.lists.pos$astyanax.markers)) {
  matrix.ast[[i]] <- lapply(gene.lists.pos$astyanax.markers, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$astyanax.markers[[i]], cell_type2 = x))
}
DI.list[[7]] <- matrix.ast

matrix.zeb.sub <- list()
for (i in 1:length(gene.lists.pos$zebrafish.markers.sub)) {
  matrix.zeb.sub[[i]] <- lapply(gene.lists.pos$zebrafish.markers.sub, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$zebrafish.markers.sub[[i]], cell_type2 = x))
}
DI.list[[8]] <- matrix.zeb.sub

matrix.ast.sub <- list()
for (i in 1:length(gene.lists.pos$astyanax.markers.sub)) {
  matrix.ast.sub[[i]] <- lapply(gene.lists.pos$astyanax.markers.sub, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$astyanax.markers.sub[[i]], cell_type2 = x))
}
DI.list[[9]] <- matrix.ast.sub

matrix.ast2zeb <- list()
for (i in 1:length(gene.lists.pos$astyanax.markers)) {
  matrix.ast2zeb[[i]] <- lapply(gene.lists.pos$zebrafish.markers, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$astyanax.markers[[i]], cell_type2 = x))
}
DI.list[[10]] <- matrix.ast2zeb

matrix.ast2zeb.sub <- list()
for (i in 1:length(gene.lists.pos$astyanax.markers.sub)) {
  matrix.ast2zeb.sub[[i]] <- lapply(gene.lists.pos$zebrafish.markers.sub, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$astyanax.markers.sub[[i]], cell_type2 = x))
}
DI.list[[11]] <- matrix.ast2zeb.sub


## Make corrected DI
# Subtract paralog numbers from species.1 and species.2 values, and add them to conserved, then calculate DI

mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")


calcCorrDriftIndex <- function(cell_type1 = cell_type1, cell_type2 = cell_type2, mart.1 = mart[[1]], mart.2 = mart[[2]]) {
  colnames(cell_type1) <- paste("cell_type1", colnames(cell_type1), sep = "_")
  colnames(cell_type2) <- paste("cell_type2", colnames(cell_type2), sep = "_")
  index <- intersect(row.names(cell_type1), row.names(cell_type2))
  df <- cbind(cell_type1[index,], cell_type2[index,])
  df$minimup_p_val <- apply(df[,grep("p_val$", colnames(df))], 1, function(x) metap::minimump(x)$p)
  df <- df[df$minimup_p_val < 0.05,]
  
  paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(row.names(df), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(row.names(df), mart.2$Gene.name)])
  paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(row.names(cell_type2), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(row.names(cell_type2), mart.1$Gene.name)])
  paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(row.names(cell_type1), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(row.names(cell_type1), mart.2$Gene.name)])
  
  paralog.1 <- paralog.1[paralog.1 %in% mart.2$Gene.name]
  paralog.2 <- paralog.2[paralog.2 %in% mart.1$Gene.name]
  
  a1 <- length(row.names(cell_type1)[row.names(cell_type1) %in% union(paralog.con, paralog.2)]) # genes.1 that are paralogs of a conserved or species 2 gene (should be subtracted from genes.1 and added to con)
  a2 <- length(row.names(cell_type2)[row.names(cell_type2) %in% union(paralog.con, paralog.1)])
  
  n.genes.1 <- nrow(cell_type1) - a1
  n.genes.2 <- nrow(cell_type2) - a2
  
  n.con <- length(row.names(df)) + a1 + a2
  
  DI.corr <- 1- sqrt( abs( (1 - (n.con / n.genes.1)) * (1 - (n.con / n.genes.2)) ) )
  DI <- 1 - sqrt( abs( (1 - (nrow(df) / nrow(cell_type1))) * (1 - (nrow(df) / nrow(cell_type2))) ) )
  return(DI.corr)
}


matrix.ast2zeb.sub.corr <- list()
for (i in 1:length(gene.lists.pos$astyanax.markers.sub)) {
  matrix.ast2zeb.sub.corr[[i]] <- lapply(gene.lists.pos$zebrafish.markers.sub, function(x) calcCorrDriftIndex(cell_type1 = gene.lists.pos$astyanax.markers.sub[[i]], cell_type2 = x))
}
names(matrix.ast2zeb.sub.corr) <- names(gene.lists.pos$astyanax.markers.sub)

DI.list[[12]] <- matrix.ast2zeb.sub.corr

## Save

names(DI.list) <- c("DI", "DI.sub", "DI.ast", "DI.ast.sub", "DI.sub.GO", "DI.zeb2zeb", "DI.ast2ast", "DI.zeb2zeb.sub", "DI.ast2ast.sub", "DI.ast2zeb", "DI.ast2zeb.sub", "DI.ast2zeb.sub.corr")

saveRDS(DI.list, file = "DI_results.rds")







## Make non row-scaled heatmap
matrix.ast2.zeb.sub <- DI.list[[11]]
names(matrix.ast2zeb.sub) <- names(matrix.ast2zeb.sub[[1]])
matrix2 <- reshape2::melt(matrix.ast2.zeb.sub)
matrix2$L1 <- unlist(lapply(names(matrix.ast2zeb.sub), function(x) rep(x, 151)))
matrix2$L1 <- factor(matrix2$L1, levels = levels(hypo.integrated@meta.data$integrated_SubclusterType))
matrix2$L2 <- factor(matrix2$L2, levels = levels(hypo.integrated@meta.data$integrated_SubclusterType))
DI.matrix <- ggplot(matrix2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_viridis(direction = -1)

## Make row-scaled heatmap
test <- lapply(matrix.ast2zeb.sub, function(x) unlist(x))
test <- lapply(test, function(x) 1-x)
test2 <- lapply(test, function(x) x/max(x))
test2 <- reshape2::melt(test2)
test2$L2 <- rep(names(matrix.ast2zeb.sub[[1]]),151)
test2$L1 <- factor(test2$L1, levels = levels(hypo.integrated@meta.data$integrated_SubclusterType))
test2$L2 <- factor(test2$L2, levels = levels(hypo.integrated@meta.data$integrated_SubclusterType))
DI.matrix.scaled <- ggplot(test2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8)) + scale_fill_viridis(direction = 1)

dev.new()
DI.matrix + DI.matrix.scaled



## Try with OG clustering

matrix.ast2zeb.og <- list()
for (i in 1:length(ast.markers[[1]])) {
  matrix.ast2zeb.og[[i]] <- lapply(zeb.markers[[1]], function(x) calcDriftIndex2(cell_type1 = ast.markers[[1]][[i]], cell_type2 = x))
}

names(matrix.ast2zeb.og) <- names(ast.markers[[1]])
matrix.ast2zeb.og.2 <- reshape2::melt(matrix.ast2zeb.og)
ggplot(matrix.ast2zeb.og.2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_viridis(direction = -1)


names(matrix.ast2zeb.sub) <- names(ast.markers[[1]])
test <- lapply(matrix.ast2zeb.sub, function(x) unlist(x))
test <- lapply(test, function(x) 1-x)
test2 <- lapply(test, function(x) x/max(x))
test2 <- reshape2::melt(test2)
test2$L2 <- rep(names(matrix.ast2zeb.sub[[1]]),151)
test2$L1 <- factor(test2$L1, levels = levels(hypo.integrated@meta.data$integrated_SubclusterType))
test2$L2 <- factor(test2$L2, levels = levels(hypo.integrated@meta.data$integrated_SubclusterType))
ggplot(test2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8)) + scale_fill_viridis(direction = 1)

test2 <- matrix.ast2zeb.og.2 %>% group_by(L1) %>% top_n(-1, value)



