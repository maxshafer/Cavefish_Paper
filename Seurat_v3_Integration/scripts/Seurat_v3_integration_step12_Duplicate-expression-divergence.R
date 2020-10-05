library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

### Functions
# Returns binary logical vector (TRUE/FALSE), match() returns first position where it matchs
# Filter the logical vector by the same logical vector

calculate.dT <- function(gene.1 = NULL, gene.2 = NULL, data = NULL) {
	gene1 <- names(unlist(lapply(data, function(x) gene.1 %in% x)))[unlist(lapply(data, function(x) gene.1 %in% x))]
	gene2 <- names(unlist(lapply(data, function(x) gene.2 %in% x)))[unlist(lapply(data, function(x) gene.2 %in% x))]
	ntu <- union(gene1, gene2)
	nti <- intersect(gene1, gene2)
	dT <- (length(ntu) - length(nti))/length(ntu)
	output <- list(ntu, nti, dT)
	output2 <- list(gene1, gene2, dT)
	return(output2)
}

no.NaNs <- function(x) {
	x[is.na(x)] <- 0
	return(x)
}

getSEQ <- function(x) {
  r <- GET(paste(server, paste("/sequence/region/Danio_rerio/", as.numeric(x$chromosome[1]),":", x[2], "..", x[3], ":1?", sep = ""), sep = ""), content_type("text/plain"))
  return(content(r))
}

### Load Seurat object and calculate average experssion per cluster (Subtype and SubclusterType)

hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")

# Load normed expression data

normed.expression <- readRDS("Normed_expression_data.rds")

str(normed.expression, max.level = 2)

norm.cluster.zeb <- normed.expression[[4]][[1]]
norm.cluster.ast <- normed.expression[[4]][[2]]

# test <- GetAssayData(hypo[[1]], assay.type = "RNA", slot = "data")
# norm.all.zeb <- data.frame(mean.exp = apply(test[,1:33000], 1, function(x) mean(x)))

norm.cluster.filtered.zeb <- lapply(norm.cluster.zeb, function(x) row.names(x)[x > 1])
norm.cluster.filtered.ast <- lapply(norm.cluster.ast, function(x) row.names(x)[x > 1])

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

nodes.zeb <- unique(mart[[1]]$Paralogue.last.common.ancestor.with.Zebrafish)
nodes.ast <- unique(mart[[2]]$Paralogue.last.common.ancestor.with.Cave.fish)

# Subset mart for each of the node levels
mart.zeb <- lapply(nodes.zeb, function(x) mart[[1]][mart[[1]]$Paralogue.last.common.ancestor.with.Zebrafish == x,])
names(mart.zeb) <- nodes.zeb
mart.ast <- lapply(nodes.ast, function(x) mart[[2]][mart[[2]]$Paralogue.last.common.ancestor.with.Cave.fish == x,])
names(mart.ast) <- nodes.ast

mart.na <- mart[[1]][!is.na(mart[[1]]$Paralogue..id..target.Zebrafish.gene.identical.to.query.gene),]
mart.id <- list()
mart.id[[1]] <- mart.na[mart.na $Paralogue..id..target.Zebrafish.gene.identical.to.query.gene < 20,]
mart.id[[2]] <- mart.na[mart.na $Paralogue..id..target.Zebrafish.gene.identical.to.query.gene > 20 & mart.na$Paralogue..id..target.Zebrafish.gene.identical.to.query.gene < 40,]
mart.id[[3]] <- mart.na[mart.na $Paralogue..id..target.Zebrafish.gene.identical.to.query.gene > 40 & mart.na$Paralogue..id..target.Zebrafish.gene.identical.to.query.gene < 60,]
mart.id[[4]] <- mart.na[mart.na $Paralogue..id..target.Zebrafish.gene.identical.to.query.gene > 60 & mart.na$Paralogue..id..target.Zebrafish.gene.identical.to.query.gene < 80,]
mart.id[[5]] <- mart.na[mart.na $Paralogue..id..target.Zebrafish.gene.identical.to.query.gene > 80,]

# Subset for only unique gene pairs (lapply doesn't work...?)

for (i in 1:length(mart.zeb)) {
	mart.zeb[[i]] <- mart.zeb[[i]][!duplicated(t(apply(mart.zeb[[i]][,c(2,4)], 1, sort))),]
}
for (i in 1:length(mart.ast)) {
	mart.ast[[i]] <- mart.ast[[i]][!duplicated(t(apply(mart.ast[[i]][,c(2,4)], 1, sort))),]
}
for (i in 1:length(mart.id)) {
	mart.id[[i]] <- mart.id[[i]][!duplicated(t(apply(mart.id[[i]][,c(2,4)], 1, sort))),]
}

# Get sorted gene pairs for each node
# These should be unique for each node

gene.pairs.zeb <- list()
for (i in 1:length(mart.zeb)) {
	gene.pairs.zeb[[i]] <- t(apply(mart.zeb[[i]][, c(2,4)], 1, sort))
}
names(gene.pairs.zeb) <- nodes.zeb

gene.pairs.ast <- list()
for (i in 1:length(mart.ast)) {
	gene.pairs.ast[[i]] <- t(apply(mart.ast[[i]][, c(2,4)], 1, sort))
}
names(gene.pairs.ast) <- nodes.ast

gene.pairs.id <- list()
for (i in 1:length(mart.id)) {
	gene.pairs.id[[i]] <- t(apply(mart.id[[i]][, c(2,4)], 1, sort))
}

# Only use those genes which are detected in my dataset (in object@raw.data)

gene.pairs.filtered.zeb <- lapply(gene.pairs.zeb, function(x) base::subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.zeb))))
gene.pairs.filtered.zeb <- lapply(gene.pairs.filtered.zeb, function(x) subset(x, x[,2] %in% unique(unlist(norm.cluster.filtered.zeb))))
# gene.pairs.filtered.zeb <- lapply(gene.pairs.filtered.zeb, function(x) subset(x, x[,1] %in% expressed.genes[[1]]))

gene.pairs.filtered.ast <- lapply(gene.pairs.ast, function(x) subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.ast))))
gene.pairs.filtered.ast <- lapply(gene.pairs.filtered.ast, function(x) subset(x, x[,2] %in% unique(unlist(norm.cluster.filtered.ast))))

gene.pairs.filtered.id <- lapply(gene.pairs.id, function(x) subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.zeb))))
gene.pairs.filtered.id <- lapply(gene.pairs.filtered.id, function(x) subset(x, x[,2] %in% unique(unlist(norm.cluster.filtered.zeb))))

# gene.pairs.filtered.zeb <- lapply(gene.pairs.zeb, function(x) subset(x, x[,1] %in% expressed.genes[[1]]))
#gene.pairs.filtered.ast <- lapply(gene.pairs.ast, function(x) subset(x, x[,1] %in% expressed.genes[[2]]))

# gene.pairs.filtered.zeb <- lapply(gene.pairs.zeb, function(x) subset(x, x[,2] %in% expressed.genes[[1]]))
# gene.pairs.filtered.ast <- lapply(gene.pairs.ast, function(x) subset(x, x[,2] %in% expressed.genes[[2]]))

# Ok for each of gene.pairs, calculate the dT

dT.zeb.out <- lapply(gene.pairs.filtered.zeb, function(y) apply(y, 1, function(x) calculate.dT(gene.1 = x[1], gene.2 = x[2], data = norm.cluster.filtered.zeb)))
dT.ast.out <- lapply(gene.pairs.filtered.ast, function(y) apply(y, 1, function(x) calculate.dT(gene.1 = x[1], gene.2 = x[2], data = norm.cluster.filtered.ast)))
dT.id.out <- lapply(gene.pairs.filtered.id, function(y) apply(y, 1, function(x) calculate.dT(gene.1 = x[1], gene.2 = x[2], data = norm.cluster.filtered.zeb)))

## Data wrangle zebrafish

dT.zeb <- lapply(dT.zeb.out, function(x) list(unlist(lapply(x, function(y) y[[3]])), lapply(x, function(y) y[[1]]), lapply(x, function(y) y[[2]])))

for (i in 1:length(dT.zeb)) {
	dT.zeb[[i]][[1]] <- no.NaNs(dT.zeb[[i]][[1]])
}

genes.dT.zeb <- list()
for (i in 2:length(dT.zeb)) {
	genes.dT.zeb[[i]] <- data.table(gene1 = gene.pairs.filtered.zeb[[i]][,1], gene2 = gene.pairs.filtered.zeb[[i]][,2], dT = as.numeric(dT.zeb[[i]][[1]]), branch = names(gene.pairs.filtered.zeb)[i])
	genes.dT.zeb[[i]]$cell.types.1 = list(dT.zeb[[i]][[2]])
	genes.dT.zeb[[i]]$cell.types.2 = list(dT.zeb[[i]][[2]])
	genes.dT.zeb[[i]]$branch <- factor(genes.dT.zeb[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Clupeocephala", "Otophysi", "Danio rerio"))
}
names(genes.dT.zeb) <- names(gene.pairs.filtered.zeb)

dTzeb <- do.call(rbind, genes.dT.zeb)

## Data wrangle astyanax

dT.ast <- lapply(dT.ast.out, function(x) list(unlist(lapply(x, function(y) y[[3]])), lapply(x, function(y) y[[1]]), lapply(x, function(y) y[[2]])))

for (i in 1:length(dT.ast)) {
	dT.ast[[i]][[1]] <- no.NaNs(dT.ast[[i]][[1]])
}

genes.dT.ast <- list()
for (i in 2:length(dT.ast)) {
	genes.dT.ast[[i]] <- data.table(gene1 = gene.pairs.filtered.ast[[i]][,1], gene2 = gene.pairs.filtered.ast[[i]][,2], dT = as.numeric(dT.ast[[i]][[1]]), branch = names(gene.pairs.filtered.ast)[i])
	genes.dT.ast[[i]]$cell.types.1 = list(dT.ast[[i]][[2]])
	genes.dT.ast[[i]]$cell.types.2 = list(dT.ast[[i]][[2]])
	genes.dT.ast[[i]]$branch <- factor(genes.dT.ast[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))
}
names(genes.dT.ast) <- names(gene.pairs.filtered.ast)

dTast <- do.call(rbind, genes.dT.ast)


dT.zeb <- lapply(dT.zeb.out, function(x) unlist(lapply(x, function(y) y[[3]])))
dT.ast <- lapply(dT.ast.out, function(x) unlist(lapply(x, function(y) y[[3]])))
dT.id <- lapply(dT.id.out, function(x) unlist(lapply(x, function(y) y[[3]])))

dT.zeb <- lapply(dT.zeb, function(x) no.NaNs(x)) # no.NaNs takes as input a vector of the dT values
dT.ast <- lapply(dT.ast, function(x) no.NaNs(x))
dT.id <- lapply(dT.id, function(x) no.NaNs(x))
lapply(dT.zeb, function(x) mean(x))
lapply(dT.ast, function(x) mean(x))
lapply(dT.id, function(x) mean(x))

genes.dT.zeb <- list()
for (i in 2:length(dT.zeb)) {
	genes.dT.zeb[[i]] <- data.frame(gene1 = gene.pairs.filtered.zeb[[i]][,1], gene2 = gene.pairs.filtered.zeb[[i]][,2], dT = as.numeric(dT.zeb[[i]]), branch = names(gene.pairs.filtered.zeb)[i])
	genes.dT.zeb[[i]]$branch <- factor(genes.dT.zeb[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Clupeocephala", "Otophysi", "Danio rerio"))
}
names(genes.dT.zeb) <- names(gene.pairs.filtered.zeb)

genes.dT.ast <- list()
for (i in 2:length(dT.ast)) {
	genes.dT.ast[[i]] <- data.frame(gene1 = gene.pairs.filtered.ast[[i]][,1], gene2 = gene.pairs.filtered.ast[[i]][,2], dT = as.numeric(dT.ast[[i]]), branch = names(gene.pairs.filtered.ast)[i])
	genes.dT.ast[[i]]$branch <- factor(genes.dT.ast[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))
}
names(genes.dT.ast) <- names(gene.pairs.filtered.ast)

genes.dT.id <- list()
for (i in 1:length(dT.id)) {
	genes.dT.id[[i]] <- data.frame(gene1 = gene.pairs.filtered.id[[i]][,1], gene2 = gene.pairs.filtered.id[[i]][,2], dT = as.numeric(dT.id[[i]]), identity = c("< 20%", "20-40%", "40-60%", "60-80%", "> 80%")[i])
}
names(genes.dT.id) <- c("< 20%", "20-40%", "40-60%", "60-80%", "> 80%")


dTzeb <- do.call(rbind, genes.dT.zeb)
dTast <- do.call(rbind, genes.dT.ast)
dTid <- do.call(rbind, genes.dT.id)

## Save dT

dT.list <- list(dTzeb, dTast, dTid)
saveRDS(dT.list, file = "dT_list.rds")



