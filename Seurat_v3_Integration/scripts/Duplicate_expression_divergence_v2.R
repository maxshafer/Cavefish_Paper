# library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(patchwork)


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
# Load Seurat v3 object with combined cluster annotations (get features from updated individual objects, integration step 1)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj") # Loads as hypo
load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_65k.Robj")


# Load normed expression data

normed.expression <- readRDS("Normed_expression_data.rds")

str(normed.expression, max.level = 2)

norm.cluster.zeb <- normed.expression[[5]][[1]]
norm.cluster.ast <- normed.expression[[5]][[2]]

# test <- GetAssayData(hypo[[1]], assay.type = "RNA", slot = "data")
# norm.all.zeb <- data.frame(mean.exp = apply(test[,1:33000], 1, function(x) mean(x)))

norm.cluster.filtered.zeb <- lapply(norm.cluster.zeb, function(x) row.names(x)[x > 1])
norm.cluster.filtered.ast <- lapply(norm.cluster.ast, function(x) row.names(x)[x > 1])

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Otophysi_Paralogs/")
mart <- list()
mart[[1]] <- read.csv("mart_export_Danio.txt", sep = "\t", head = TRUE)
mart[[2]] <- read.csv("~/Downloads/mart_export_Astyanax_102.txt", sep = "\t", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")
# Node levels for phylo tree for the origin of duplicated gene pairs
nodes.zeb <- levels(mart[[1]]$Paralogue.last.common.ancestor.with.Zebrafish)
nodes.ast <- levels(mart[[2]]$Paralogue.last.common.ancestor.with.Cave.fish)
            

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
# expressed.genes <- list()
# expressed.genes[[1]] <- rbindlist(norm.cluster.filtered.zeb)
# expressed.genes[[2]] <- rbindlist(norm.cluster.filtered.ast)
# expressed.genes <- lapply(expressed.genes, function(x) unique(x$gene))
# gene.pairs.filtered.zeb <- lapply(gene.pairs.zeb, function(x) subset(x, x[,1] %in% expressed.genes[[1]] | x[,2] %in% expressed.genes[[1]]))
# gene.pairs.filtered.ast <- lapply(gene.pairs.ast, function(x) subset(x, x[,1] %in% expressed.genes[[2]] | x[,2] %in% expressed.genes[[2]]))

gene.pairs.filtered.zeb <- lapply(gene.pairs.zeb, function(x) subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.zeb))))
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
	genes.dT.ast[[i]]$branch <- factor(genes.dT.ast[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Clupeocephala", "Otophysi", "Danio rerio"))
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

## Change factor levels for plotting colors



plot.zeb <- ggplot(dTzeb, aes(x = dT, color = branch)) + stat_ecdf(size = 1) + theme(axis.text.x = element_blank(), axis.title = element_blank())
plot.ast <- ggplot(dTast, aes(x = dT, color = branch)) + stat_ecdf(size = 1) + ylab("Empirical cumulative distribution") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
plot.id <- ggplot(dTid, aes(x = dT, color = identity)) + stat_ecdf(size = 1) + ylab("Empirical cumulative distribution") + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10))

plot.zeb / plot.ast / plot.id + plot_annotation(title = "Duplicate Divergence by gene age or % identity")


## Compare gene pairs which are common to zebrafish and cavefish
## Compare dTs for cavefish and zebrafish, and correlate binarized gene expression patterns
library(plyr)

dTzeb.2 <- match_df(dTzeb, dTast, on = c("gene1", "gene2"))
dTast.2 <- match_df(dTast, dTzeb, on = c("gene1", "gene2"))
dT.combined <- merge(dTast.2, dTzeb.2, by = c("gene1", "gene2"))

ggplot(dT.combined, aes(x = dT.x, y = dT.y, color = branch.x, label = paste(gene1, gene2, sep = "-"))) + geom_point(shape = 21) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + facet_wrap(~branch.x, scales = "free") + guides(color = F)


unique(dT.combined$branch.x)
[1] Opisthokonta  Bilateria     <NA>          Euteleostomi  Vertebrata    Clupeocephala Chordata      Neopterygii   Otophysi 


dTcombined <- lapply(unique(dT.combined$branch.x), function(x) dT.combined[dT.combined$branch.x == x,])
names(dTcombined) <- unique(dT.combined$branch.x)


test <- lapply(dTcombined, function(y) apply(y, 1, function(x) cor(data.table(para1 = as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[5]))), para2 = as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[6])))), data.table(para1 = as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[9]))), para2 = as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[10])))))))

test <- lapply(dTcombined, function(y) apply(y, 1, function(x) cor(c(as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[5]))), as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[6])))), c(as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[9]))), as.numeric(integrated_SubclusterTypes %in% unlist(as.list(x[10])))))))

lapply(test, function(x) mean(x))

test2 <- data.frame(branch = c(rep("Opisthokonta", length(test[[1]])), rep("Bilateria", length(test[[2]])), rep("none", length(test[[3]])), rep("Euteleostomi", length(test[[4]])), rep("Vertebrata", length(test[[5]])), rep("Clupeocephala", length(test[[6]])), rep("Chordata", length(test[[7]])), rep("Neopterygii", length(test[[8]])), rep("Otophysi", length(test[[9]]))), value = unlist(test))

test2$branch <- factor(test2$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))

anno_df <- compare_means(value ~ branch, data = test2)

ggplot(test2, aes(x = branch, y = value)) + geom_jitter(color = "grey65", shape = 21) + geom_boxplot(fill = "transparent", outlier.color = NA) + theme(axis.text = element_text(size = 8), axis.title = element_text(size = 10), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + stat_compare_means()









Idents(hypo.integrated) <- "species.2"

DimPlot(hypo.integrated, reduction = "tsne", group.by = "integrated_Subtype", label = T) + NoLegend()
FeaturePlot(hypo.integrated, features = c("nebl", "zic2a"), reduction = "tsne", min.cutoff = "q9")

DotPlot(subset(hypo.integrated, idents = "zebrafish"), features = c("ckba", "ckbb"), group.by = "integrated_Subtype")
DotPlot(subset(hypo.integrated, idents = "astyanax"), features = c("tmem107", "tmem107l"), group.by = "integrated_SubclusterType") + coord_flip()



DotPlot(subset(hypo.integrated, idents = "zebrafish"), features = c("vmo1a", "vmo1b"), group.by = "integrated_SubclusterType") + coord_flip()
DotPlot(subset(hypo.integrated, idents = "astyanax"), features = c("stmn2a", "stmn2b"), group.by = "integrated_SubclusterType") + coord_flip()

DotPlot(hypo.integrated, features = c("nebl"), group.by = "integrated_Subtype", split.by = "species.2", cols = c("blue", "red"))


DotPlot(subset(hypo.integrated, idents = "astyanax"), features = c("six6b"), group.by = "integrated_Subtype", split.by = "species")

FeaturePlot(hypo.ast, features = c("qkia"), reduction = "FItSNE", min.cutoff = "q9")
DotPlot(hypo.zeb, features = c("gria3a", "grin1a"), group.by = "Subtype")

