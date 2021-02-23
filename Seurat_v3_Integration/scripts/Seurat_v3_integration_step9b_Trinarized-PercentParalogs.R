library(dplyr)
library(fdrtool)

# Set hyperparameters from arguments or in file
args <- commandArgs(trailingOnly = TRUE)

print(args)

a <- as.numeric(args[[1]])
b <- as.numeric(args[[2]])
f <- as.numeric(args[[3]])

# a = 1.5
# b = 2
# f = 0.35

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

## Take sets of conserved, and species specific marker genes (start with Subtypes)
## For each species, take each gene (apply), and query mart database for paralogs
## Ask if any of the paralogs are in species.2 list
## Append species.1 list with paralog expressed (if more then 1??), plus LCA


### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

## Load trinarized gene lists

trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))


### Association between whether a gene is species-expressed in a cell type, and it being a paralog of gene expressed in either species, or conserved between species.

## Identify all paralogs of all expressed genes union(paralogs of conserved, paralogs of species.1, paraglos of species.2)
## Ask what % of species specific expressed genes are in the above list

calcParalog <- function(species.1 = species.1, species.2 = species.2, mart.1 = mart[[1]], mart.2 = mart[[2]], ngenes.1 = 32191, ngenes.2 = 25271, i = i) {
  names <- Reduce(intersect, list(names(species.1), names(species.2)))
  
  conserved <- lapply(names, function(x) intersect(names(species.1[[x]]), names(species.2[[x]])))
  species.1 <- lapply(names, function(x) names(species.1[[x]]))
  species.2 <- lapply(names, function(x) names(species.2[[x]]))
  
  genes.1 <- setdiff(species.1[[i]], conserved[[i]])
  genes.2 <- setdiff(species.2[[i]], conserved[[i]])
  
  paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(conserved[[i]], mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(conserved[[i]], mart.2$Gene.name)])
  paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart.1$Gene.name)])
  paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart.2$Gene.name)])
  paralog.union <- union(paralog.con, union(paralog.1, paralog.2))
  
  paralog.1 <- paralog.1[paralog.1 %in% mart.2$Gene.name]
  paralog.2 <- paralog.2[paralog.2 %in% mart.1$Gene.name]
  
  a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.2)])
  b1 <- length(genes.1) - a1
  c1 <- length(union(paralog.con, paralog.2)) - a1
  d1 <- ngenes.1 - b1
  
  a2 <- length(genes.2[genes.2 %in% union(paralog.con, paralog.1)])
  b2 <- length(genes.2) - a2
  c2 <- length(union(paralog.con, paralog.1)) - a2
  d2 <- ngenes.2 - b2
  
  a3 <- length(row.names(conserved[[i]]))
  b3 <- length(genes.1)
  c3 <- length(genes.2)
  d3 <- ngenes.1 - sum(b3, c3)
  
  vec <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2, a3=a3, b3=b3, c3=c3, d3=d3)
  return(vec)
}


paralog.numbers <- lapply(seq_along(trinarized.exp[[1]]), function(x) calcParalog(species.1 = trinarized.exp[[1]], species.2 = trinarized.exp[[2]], i = x))

fisher.results.1 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(trinarized.exp[[1]])
fisher.results.2 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(trinarized.exp[[1]])
fisher.results.3 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2)))
fisher.results.3 <- do.call(rbind, fisher.results.3)
row.names(fisher.results.3) <- names(trinarized.exp[[1]])

clusters <- data.frame(do.call(rbind, paralog.numbers))
row.names(clusters) <- names(trinarized.exp[[1]])
clusters$percent.para.1 <- unlist(apply(clusters, 1, function(x) x[1]/sum(x[1], x[3])))*100
clusters$percent.para.2 <- unlist(apply(clusters, 1, function(x) x[5]/sum(x[5], x[7])))*100
# clusters$fisher.pval <- as.numeric(fisher.results[,1])
clusters$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
clusters$fisher.pvalue.1 <- as.numeric(fisher.results.1[,1])
clusters$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
clusters$fisher.pvalue.2 <- as.numeric(fisher.results.2[,1])
clusters$Clusters <- row.names(clusters)
clusters$similarityindex <- SI.list[[1]]$values

clusters$Clusters <- factor(clusters$Clusters, levels = levels(hypo.integrated@meta.data$integrated_Cluster))



#### For subclusters!

paralog.numbers.sub <- lapply(seq_along(trinarized.exp[[5]]), function(x) calcParalog(species.1 = trinarized.exp[[5]], species.2 = trinarized.exp[[6]], i = x))

fisher.results.1 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(trinarized.exp[[5]])
fisher.results.2 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(trinarized.exp[[5]])

subclusters <- data.frame(do.call(rbind, paralog.numbers.sub))
row.names(subclusters) <- names(trinarized.exp[[5]])
subclusters$percent.para.1 <- unlist(apply(subclusters, 1, function(x) x[1]/sum(x[1], x[3])))*100
subclusters$percent.para.2 <- unlist(apply(subclusters, 1, function(x) x[5]/sum(x[5], x[7])))*100
# subclusters$fisher.pval <- as.numeric(fisher.results[,1])
subclusters$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
subclusters$fisher.pvalue.1 <- as.numeric(fisher.results.1[,1])
subclusters$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
subclusters$fisher.pvalue.2 <- as.numeric(fisher.results.2[,1])
subclusters$Subclusters <- row.names(subclusters)
subclusters$similarityindex <- SI.list[[2]]$values
subclusters$Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(subclusters$Subclusters, hypo.integrated@meta.data$integrated_Subcluster)]

subclusters$Cluster <- factor(subclusters$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

## Save results

para.results <- list(subtypes, subclusters)

saveRDS(para.results, file = paste("Paralog-results-trinarized_a",a,"_b",b, "_f",f,".rds", sep = ""))


