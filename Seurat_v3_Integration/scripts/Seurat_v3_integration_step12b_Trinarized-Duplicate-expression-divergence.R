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

hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")

zeb.genes <- rownames(GetAssayData(hypo.zeb))
ast.genes <- rownames(GetAssayData(hypo.ast))

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

# Only use those genes which are detected in my dataset (in trinarized genes per cluster)

## Load trinarized gene lists

trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# Make sure to include both the subclusters, and the species-specific subclusters
norm.cluster.filtered.zeb <- lapply(c(trinarized.exp$subcluster.zebrafish, trinarized.exp$specific.zebrafish), function(x) names(x))
norm.cluster.filtered.ast <- lapply(c(trinarized.exp$subcluster.astyanax, trinarized.exp$specific.astyanax), function(x) names(x))


gene.pairs.filtered.zeb <- lapply(gene.pairs.zeb, function(x) base::subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.zeb)) | x[,2] %in% unique(unlist(norm.cluster.filtered.zeb))))

gene.pairs.filtered.ast <- lapply(gene.pairs.ast, function(x) base::subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.ast)) | x[,2] %in% unique(unlist(norm.cluster.filtered.ast))))

gene.pairs.filtered.id <- lapply(gene.pairs.id, function(x) base::subset(x, x[,1] %in% unique(unlist(norm.cluster.filtered.zeb)) | x[,2] %in% unique(unlist(norm.cluster.filtered.zeb))))

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
saveRDS(dT.list, file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))



## Find the ancestral cell types and use them to compute ancestral dT and NF scores
## ancestral == intersect(cell_types(GeneXa | GeneXb), cell_types(GeneXa  | GeneXb))
## NF score == # of new cell types / total cell types

## This makes a structured list of lists (2x (2, 2, 1)) of cell types for each gene, in each species + ancestral expression pattern
## Use a separate function to calculate the dT and NF values (all permutations)

calcCellTypes <- function(gene.pairs = gene.pairs, species.1.genes = zeb.genes, species.2.genes = ast.genes, data = list(norm.cluster.filtered.zeb, norm.cluster.filtered.ast)) {
  
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

calcAncestraldT <- function(cell.types = cell.types.list.element, sp1 = length(norm.cluster.filtered.zeb), sp2 = length(norm.cluster.filtered.ast)) {
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
colnames(dT.combined2) <- c("gene1", "gene2", "dT.zeb", "branch.zeb", "dT.ast", "branch.ast")

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
colnames(dT.NF) <- c("gene1","gene2","dT.zeb","branch.zeb","dT.ast","branch.ast","dT.zeb.ancestral","dT.ast.ancestral","NF.gene1.zeb", "NF.gene2.zeb", "NF.gene1.ast", "NF.gene2.ast")

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
# gene1  gene2 dT.zeb branch.zeb dT.ast   branch.ast dT.zeb.ancestral dT.ast.ancestral NF.gene1.zeb NF.gene2.zeb NF.gene1.ast NF.gene2.ast       branch
# 1        aadac     NA       <NA>      1     Chordata               NA               NA           NA   0.00000000           NA   0.01986755     Chordata
# 2       abcb6a     NA       <NA>      1    Bilateria               NA               NA           NA   0.00000000           NA   0.01986755    Bilateria
# 3        abcc3     NA       <NA>      1    Bilateria               NA               NA           NA   0.00000000           NA   0.03311258    Bilateria
# 4        abcc9     NA       <NA>      1 Euteleostomi               NA               NA           NA   0.00000000           NA   0.01324503 Euteleostomi
# 5        abcc9     NA       <NA>      1   Vertebrata               NA               NA           NA   0.00000000           NA   0.01324503   Vertebrata
# 6       abhd12      1 Vertebrata      1   Vertebrata               NA               NA           NA   0.01986755           NA   0.04635762   Vertebrata


## Save dT list

dT.list[[4]] <- dT.NF
saveRDS(dT.list, file = paste("dT_list-trinarized_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

