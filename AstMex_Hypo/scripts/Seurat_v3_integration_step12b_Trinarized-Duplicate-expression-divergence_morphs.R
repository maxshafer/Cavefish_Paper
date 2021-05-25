library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(data.table)
library(patchwork)

# Set hyperparameters from arguments or in file
# args <- commandArgs(trailingOnly = TRUE)
# 
# print(args)
# 
# a <- as.numeric(args[[1]])
# b <- as.numeric(args[[2]])
# f <- as.numeric(args[[3]])

a = 1.5
b = 2
f = 0.2

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

### Functions
# Returns binary logical vector (TRUE/FALSE), match() returns first position where it matchs
# Filter the logical vector by the same logical vector

calculate.dT <- function(gene.1 = NULL, gene.2 = NULL, data = NULL) {
	gene1 <- names(unlist(lapply(data, function(x) gene.1 %in% x)))[unlist(lapply(data, function(x) gene.1 %in% x))] # This is the cell types where gene.1 appears
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

hypo.cave <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_vR.rds")

ast.genes <- rownames(GetAssayData(hypo.ast))

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

nodes.ast <- unique(mart[[2]]$Paralogue.last.common.ancestor.with.Cave.fish)

# Subset mart for each of the node levels
mart.ast <- lapply(nodes.ast, function(x) mart[[2]][mart[[2]]$Paralogue.last.common.ancestor.with.Cave.fish == x,])
names(mart.ast) <- nodes.ast

# Subset for only unique gene pairs (lapply doesn't work...?)

for (i in 1:length(mart.ast)) {
	mart.ast[[i]] <- mart.ast[[i]][!duplicated(t(apply(mart.ast[[i]][,c(2,4)], 1, sort))),]
}

# Get sorted gene pairs for each node
# These should be unique for each node

gene.pairs.ast <- list()
for (i in 1:length(mart.ast)) {
	gene.pairs.ast[[i]] <- t(apply(mart.ast[[i]][, c(2,4)], 1, sort))
}
names(gene.pairs.ast) <- nodes.ast

# Only use those genes which are detected in my dataset (in trinarized genes per cluster)

## Load trinarized gene lists

trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# Make sure to include both the subclusters, and the species-specific subclusters
norm.cluster.filtered.surface <- lapply(trinarized.exp$subcluster.surface, function(x) names(x))
norm.cluster.filtered.cave <- lapply(trinarized.exp$subcluster.cave, function(x) names(x))

gene.pairs.filtered.surface <- lapply(gene.pairs.ast, function(x) base::subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.surface)) & x[,2] %in% unique(unlist(norm.cluster.filtered.surface))))
gene.pairs.filtered.cave <- lapply(gene.pairs.ast, function(x) base::subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.cave)) & x[,2] %in% unique(unlist(norm.cluster.filtered.cave))))

# Ok for each of gene.pairs, calculate the dT

dT.surface.out <- lapply(gene.pairs.filtered.surface, function(y) apply(y, 1, function(x) calculate.dT(gene.1 = x[1], gene.2 = x[2], data = norm.cluster.filtered.surface)))
dT.cave.out <- lapply(gene.pairs.filtered.cave, function(y) apply(y, 1, function(x) calculate.dT(gene.1 = x[1], gene.2 = x[2], data = norm.cluster.filtered.cave)))

## Data wrangle zebrafish

dT.surface <- lapply(dT.surface.out, function(x) list(unlist(lapply(x, function(y) y[[3]])), lapply(x, function(y) y[[1]]), lapply(x, function(y) y[[2]])))

for (i in 1:length(dT.surface)) {
	dT.surface[[i]][[1]] <- no.NaNs(dT.surface[[i]][[1]])
}

genes.dT.surface <- list()
for (i in 2:length(dT.surface)) {
	genes.dT.surface[[i]] <- data.table(gene1 = gene.pairs.filtered.surface[[i]][,1], gene2 = gene.pairs.filtered.surface[[i]][,2], dT = as.numeric(dT.surface[[i]][[1]]), branch = names(gene.pairs.filtered.surface)[i])
	genes.dT.surface[[i]]$cell.types.1 = list(dT.surface[[i]][[2]])
	genes.dT.surface[[i]]$cell.types.2 = list(dT.surface[[i]][[2]])
	genes.dT.surface[[i]]$branch <- factor(genes.dT.surface[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Clupeocephala", "Otophysi", "Danio rerio"))
}
names(genes.dT.surface) <- names(gene.pairs.filtered.surface)

dTsurface <- do.call(rbind, genes.dT.surface)

## Data wrangle astyanax

dT.cave <- lapply(dT.cave.out, function(x) list(unlist(lapply(x, function(y) y[[3]])), lapply(x, function(y) y[[1]]), lapply(x, function(y) y[[2]])))

for (i in 1:length(dT.cave)) {
	dT.cave[[i]][[1]] <- no.NaNs(dT.cave[[i]][[1]])
}

genes.dT.cave <- list()
for (i in 2:length(dT.cave)) {
	genes.dT.cave[[i]] <- data.table(gene1 = gene.pairs.filtered.cave[[i]][,1], gene2 = gene.pairs.filtered.cave[[i]][,2], dT = as.numeric(dT.cave[[i]][[1]]), branch = names(gene.pairs.filtered.cave)[i])
	genes.dT.cave[[i]]$cell.types.1 = list(dT.cave[[i]][[2]])
	genes.dT.cave[[i]]$cell.types.2 = list(dT.cave[[i]][[2]])
	genes.dT.cave[[i]]$branch <- factor(genes.dT.cave[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))
}
names(genes.dT.cave) <- names(gene.pairs.filtered.cave)

dTcave <- do.call(rbind, genes.dT.cave)


dT.surface <- lapply(dT.surface.out, function(x) unlist(lapply(x, function(y) y[[3]])))
dT.cave <- lapply(dT.cave.out, function(x) unlist(lapply(x, function(y) y[[3]])))

dT.surface <- lapply(dT.surface, function(x) no.NaNs(x)) # no.NaNs takes as input a vector of the dT values
dT.cave <- lapply(dT.cave, function(x) no.NaNs(x))
lapply(dT.surface, function(x) mean(x))
lapply(dT.cave, function(x) mean(x))

genes.dT.surface <- list()
for (i in 2:length(dT.surface)) {
	genes.dT.surface[[i]] <- data.frame(gene1 = gene.pairs.filtered.surface[[i]][,1], gene2 = gene.pairs.filtered.surface[[i]][,2], dT = as.numeric(dT.surface[[i]]), branch = names(gene.pairs.filtered.surface)[i])
	genes.dT.surface[[i]]$branch <- factor(genes.dT.surface[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Clupeocephala", "Otophysi", "Danio rerio"))
}
names(genes.dT.surface) <- names(gene.pairs.filtered.surface)

genes.dT.cave <- list()
for (i in 2:length(dT.cave)) {
	genes.dT.cave[[i]] <- data.frame(gene1 = gene.pairs.filtered.cave[[i]][,1], gene2 = gene.pairs.filtered.cave[[i]][,2], dT = as.numeric(dT.cave[[i]]), branch = names(gene.pairs.filtered.cave)[i])
	genes.dT.cave[[i]]$branch <- factor(genes.dT.cave[[i]]$branch, levels = c("Opisthokonta", "Bilateria", "Chordata", "Vertebrata", "Euteleostomi", "Neopterygii", "Osteoglossocephalai", "Clupeocephala", "Otophysi", "Characiphysae", "Characoidei", "Astyanax mexicanus"))
}
names(genes.dT.cave) <- names(gene.pairs.filtered.cave)


dTsurface <- do.call(rbind, genes.dT.surface)
dTcave <- do.call(rbind, genes.dT.cave)

## Save dT

dT.list <- list(dTsurface, dTcave)

saveRDS(dT.list, file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

dTzeb.2 <- plyr::match_df(dT.list[[1]], dT.list[[2]], on = c("gene1", "gene2"))
dTast.2 <- plyr::match_df(dT.list[[2]], dT.list[[1]], on = c("gene1", "gene2"))
dT.combined <- merge(dTast.2, dTzeb.2, by = c("gene1", "gene2"))

# x is ast, y is zeb
all.genes <- ggplot(dT.combined, aes(x = dT.x, y = dT.y, color = branch.x, label = paste(gene1, gene2, sep = "-"))) + geom_point(colour = "black", size = 0.5) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + theme_classic() + ylab("dT zebrafish") + xlab("dT Mexican tetra") + ggtitle("All gene pairs", subtitle = "Pearson correlation = 0.54")
otophysi <- ggplot(dT.combined[dT.combined$branch.x == "Otophysi" | dT.combined$branch.y == "Otophysi",], aes(x = dT.x, y = dT.y)) + geom_point(colour = "#FF61C9", size = 0.5) + theme(axis.text = element_text(size = 8), strip.text = element_text(size = 8)) + theme_classic() + ylab("dT zebrafish") + xlab("dT Mexican tetra") + ggtitle("LCA Otophysi gene pairs", subtitle = "Pearson correlation = 0.45")


## Find the ancestral cell types and use them to compute ancestral dT and NF scores
## ancestral == intersect(cell_types(GeneXa | GeneXb), cell_types(GeneXa  | GeneXb))
## NF score == # of new cell types / total cell types

## This makes a structured list of lists (2x (2, 2, 1)) of cell types for each gene, in each species + ancestral expression pattern
## Use a separate function to calculate the dT and NF values (all permutations)

calcCellTypes <- function(gene.pairs = gene.pairs, species.1.genes = zeb.genes, species.2.genes = ast.genes, data = list(norm.cluster.filtered.surface, norm.cluster.filtered.cave)) {
  
  species.genes <- list(species.1.genes, species.2.genes)
  # For each pair of genes (supplied by lapply), find the cell types that express each in both species (data)
  cell.types <- list()
  for (i in 1:length(data)) {
    cell.types[[i]] <- list()
    for (j in 1:length(gene.pairs)) {
      if (gene.pairs[[j]] %in% species.genes[[i]]) {
        cell.types[[i]][[j]] <- names(unlist(lapply(data[[i]], function(x) gene.pairs[[j]] %in% x)))[unlist(lapply(data[[i]], function(x) gene.pairs[[j]] %in% x))]
      } else {
        cell.types[[i]][[j]] <- "no_gene"
      }
    }
  }
  
  ancestral <- intersect(Reduce(union, cell.types[[1]][cell.types[[1]] != "no_gene"]),Reduce(union, cell.types[[2]][cell.types[[2]] != "no_gene"]))
  if (is.null(ancestral)) {
    ancestral <- "both genes are species-specific"
  }
  
  cell.types[[3]] <- ancestral
  names(cell.types) <- c("zebrafish", "mexican_tetra", "ancestral")
  names(cell.types[[1]]) <- c(gene.pairs[[1]], gene.pairs[[2]])
  names(cell.types[[2]]) <- c(gene.pairs[[1]], gene.pairs[[2]])
  return(cell.types)
}

calcAncestraldT <- function(cell.types = cell.types.list.element, sp1 = length(norm.cluster.filtered.surface), sp2 = length(norm.cluster.filtered.cave)) {
  # Calculate scores as long as it isn't "no_gene", if so, put NA
  # cell.types is a list length 3, each of the first two (for species) are lists of length 2 for each gene, 3 is ancestral
  sp <- c(sp1,sp2)
  # Calculate NF score
  nf.scores <- lapply(seq_along(cell.types[1:2]), function(x) lapply(cell.types[[x]], function(y) {
    if("both genes are species-specific" %in% cell.types[[3]]) {
      return(NA)
    } else {
      if("no_gene" %in% y){
        return(NA)
      }
    }
    length(y[!(y %in% cell.types[[3]])])/sp[[x]]
  }))
  names(nf.scores) <- c("NF_zebrafish", "NF_mexican_tetra")
  
  anc.cell.types <- lapply(cell.types[1:2], function(x) lapply(x, function(y) y[y %in% cell.types[[3]]]))
  
  ntu <- lapply(anc.cell.types, function(x) length(union(x[[1]], x[[2]])))
  nti <- lapply(anc.cell.types, function(x) length(intersect(x[[1]], x[[2]])))
  
  dT <- lapply(seq_along(ntu), function(x) {
    if (is.na((ntu[[x]] - nti[[x]])/ntu[[x]])) {
      return(NA)
    } else {
      (ntu[[x]] - nti[[x]])/ntu[[x]]
    }
  })
  names(dT) <- c("dT_zebrafish", "dT_mexican_tetra")
  
  dT[unlist(lapply(cell.types[1:2], function(x) "no_gene" %in% x))] <- NA
  dT[unlist(lapply(cell.types[1:2], function(x) "no_comparison" %in% x))] <- NA
  
  output2 <- c(unlist(dT), unlist(nf.scores))
  return(output2)
}

# Set variables and load dT list

# a = 1.5
# b = 2
# f = 0.35

dT.list <- readRDS(file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

dT.combined2 <- merge(dT.list[[1]], dT.list[[2]], by = c("gene1", "gene2"), all = T)
colnames(dT.combined2) <- c("gene1", "gene2", "dT.surface", "branch.surface", "dT.cave", "branch.cave")

## Now run calcCellTypes and calcAncestraldT on the dT.combined2

# Calculate CellTypes
CellTypes <- apply(dT.combined2, 1, function(x) {
  return(tryCatch(calcCellTypes(gene.pairs = x[1:2]), error=function(e) NULL))
})
names(CellTypes) <- paste(dT.combined2[,1], dT.combined2[,2], sep = "_")

dT.NF <- lapply(CellTypes, function(x) calcAncestraldT(cell.types = x))
dT.NF <- Reduce(rbind, dT.NF)
rownames(dT.NF) <- rownames(dT.combined2)
dT.NF <- cbind(dT.combined2, dT.NF)
colnames(dT.NF) <- c("gene1","gene2","dT.surface","branch.surface","dT.cave","branch.cave","dT.surface.ancestral","dT.cave.ancestral","NF.gene1.surface", "NF.gene2.surface", "NF.gene1.cave", "NF.gene2.cave")

branch.levels <- levels(dT.list[[2]]$branch)
branch.levels <- c(branch.levels, "Danio rerio")

dT.NF$branch <- factor(unlist(apply(dT.NF, 1, function(x) {
  name <- unique(x[c(4,6)][!(is.na(x[c(4,6)]))])
  if (length(name) < 1) {
    name <- "empty"
  }
  if (length(name) > 1) {
    name <- paste(name[[1]], name[[2]], sep = "_")
  }
  return(name)
}
)), levels = branch.levels)

# dT.NF is now a df with rows for each gene pair (where at least 1 member of the pair was *expressed* in at least 1 subcluster in either species)
# and columns for original dT scores, ancestral dT scores, and NF scores for each gene (or NAs for missing values where applicable)
# Inclusion of pairs where only 1 gene is detected does inflate the dT = 1 values (since these are by definition diverged in the hypothalamus) 

head(dT.NF)
# gene1  gene2 dT.surface branch.surface dT.cave   branch.cave dT.surface.ancestral dT.cave.ancestral NF.gene1.surface NF.gene2.surface NF.gene1.cave NF.gene2.cave       branch
# 1        aadac     NA       <NA>      1     Chordata               NA               NA           NA   0.00000000           NA   0.01986755     Chordata
# 2       abcb6a     NA       <NA>      1    Bilateria               NA               NA           NA   0.00000000           NA   0.01986755    Bilateria
# 3        abcc3     NA       <NA>      1    Bilateria               NA               NA           NA   0.00000000           NA   0.03311258    Bilateria
# 4        abcc9     NA       <NA>      1 Euteleostomi               NA               NA           NA   0.00000000           NA   0.01324503 Euteleostomi
# 5        abcc9     NA       <NA>      1   Vertebrata               NA               NA           NA   0.00000000           NA   0.01324503   Vertebrata
# 6       abhd12      1 Vertebrata      1   Vertebrata               NA               NA           NA   0.01986755           NA   0.04635762   Vertebrata


## Save dT list

dT.list[[4]] <- dT.NF
saveRDS(dT.list, file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

