library(dplyr)
library(fdrtool)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

## Take sets of conserved, and species specific marker genes (start with Clusters)
## For each species, take each gene (apply), and query mart database for paralogs
## Ask if any of the paralogs are in species.2 list
## Append species.1 list with paralog expressed (if more then 1??), plus LCA


### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

## Load marker gene lists

a = 1.5
b = 2
f = 0.35

gene.lists.pos <- readRDS(paste("drift_gene_lists_pos_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))
SI.list <- readRDS(paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))


### Association between whether a gene is a specieis specific marker, and it being a paralog of a marker gene in either species, or conserved between species.

## Identify all paralogs of all marker genes union(paralogs of conserved, paralogs of species.1, paraglos of species.2)
## Ask what % of species specific marker genes are in the above list

calcParalog <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, mart.1 = mart[[1]], mart.2 = mart[[2]], ngenes.1 = 32191, ngenes.2 = 25271, i = i) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	
	genes.1 <- setdiff(row.names(species.1[[i]]), row.names(conserved[[i]]))
	genes.2 <- setdiff(row.names(species.2[[i]]), row.names(conserved[[i]]))
	
	paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(row.names(conserved[[i]]), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(row.names(conserved[[i]]), mart.2$Gene.name)])
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


paralog.numbers <- lapply(seq_along(gene.lists.pos[[1]]), function(x) calcParalog(conserved = gene.lists.pos[["cluster.conserved"]], species.1 = gene.lists.pos[["cluster.zebrafish"]], species.2 = gene.lists.pos[["cluster.astyanax"]], i = x))

fisher.results.1 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[[1]])
fisher.results.2 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(gene.lists.pos[[1]])
fisher.results.3 <- lapply(paralog.numbers, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2)))
fisher.results.3 <- do.call(rbind, fisher.results.3)
row.names(fisher.results.3) <- names(gene.lists.pos[[1]])

clusters <- data.frame(do.call(rbind, paralog.numbers))
row.names(clusters) <- names(gene.lists.pos[[1]])
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

paralog.numbers.sub <- lapply(seq_along(1:151), function(x) calcParalog(conserved = gene.lists.pos[["subcluster.conserved"]], species.1 = gene.lists.pos[["subcluster.zebrafish"]], species.2 = gene.lists.pos[["subcluster.astyanax"]], i = x))

fisher.results.1 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[["subcluster.conserved"]])
fisher.results.2 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(gene.lists.pos[[7]])

subclusters <- data.frame(do.call(rbind, paralog.numbers.sub))
row.names(subclusters) <- names(gene.lists.pos[[7]])
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

#### Do the same but for cave vs surface

paralog.numbers.ast <- lapply(seq_along(gene.lists.pos[[1]]), function(x) calcParalog(conserved = gene.lists.pos[["cluster.conserved.ast"]], species.1 = gene.lists.pos[["cluster.surface"]], species.2 = gene.lists.pos[["cluster.cave"]], mart.1 = mart[[2]], mart.2 = mart[[2]], i = x))

fisher.results.1 <- lapply(paralog.numbers.ast, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[[1]])
fisher.results.2 <- lapply(paralog.numbers.ast, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(gene.lists.pos[[1]])
fisher.results.3 <- lapply(paralog.numbers.ast, function(x) fisher.test(matrix(data = x[9:12], ncol = 2, nrow = 2)))
fisher.results.3 <- do.call(rbind, fisher.results.3)
row.names(fisher.results.3) <- names(gene.lists.pos[[1]])

clusters.ast <- data.frame(do.call(rbind, paralog.numbers.ast))
row.names(clusters.ast) <- names(gene.lists.pos[[1]])
clusters.ast$percent.para.1 <- unlist(apply(clusters.ast, 1, function(x) x[1]/sum(x[1], x[3])))*100
clusters.ast$percent.para.2 <- unlist(apply(clusters.ast, 1, function(x) x[5]/sum(x[5], x[7])))*100
# clusters.ast$fisher.pval <- as.numeric(fisher.results[,1])
clusters.ast$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
clusters.ast$fisher.pvalue.1 <- as.numeric(fisher.results.1[,1])
clusters.ast$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
clusters.ast$fisher.pvalue.2 <- as.numeric(fisher.results.2[,1])
clusters.ast$Clusters <- row.names(clusters.ast)
clusters.ast$similarityindex <- SI.list[[1]]$values

clusters.ast$Clusters <- factor(clusters.ast$Clusters, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

#### For subclusters!

paralog.numbers.sub.ast <- lapply(seq_along(gene.lists.pos[[7]]), function(x) calcParalog(conserved = gene.lists.pos[["subcluster.conserved.ast"]], species.1 = gene.lists.pos[["subcluster.surface"]], species.2 = gene.lists.pos[["subcluster.cave"]], mart.1 = mart[[2]], mart.2 = mart[[2]], i = x))

fisher.results.1 <- lapply(paralog.numbers.sub.ast, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[[7]])
fisher.results.2 <- lapply(paralog.numbers.sub.ast, function(x) fisher.test(matrix(data = x[5:8], ncol = 2, nrow = 2)))
fisher.results.2 <- do.call(rbind, fisher.results.2)
row.names(fisher.results.2) <- names(gene.lists.pos[[7]])

subclusters.ast <- data.frame(do.call(rbind, paralog.numbers.sub.ast))
row.names(subclusters.ast) <- names(gene.lists.pos[[7]])
subclusters.ast$percent.para.1 <- unlist(apply(subclusters.ast, 1, function(x) x[1]/sum(x[1], x[3])))*100
subclusters.ast$percent.para.2 <- unlist(apply(subclusters.ast, 1, function(x) x[5]/sum(x[5], x[7])))*100
# subclusters.ast$fisher.pval <- as.numeric(fisher.results[,1])
subclusters.ast$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
subclusters.ast$fisher.pvalue.1 <- as.numeric(fisher.results.1[,1])
subclusters.ast$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
subclusters.ast$fisher.pvalue.2 <- as.numeric(fisher.results.2[,1])
subclusters.ast$Subclusters <- row.names(subclusters.ast)
subclusters.ast$similarityindex <- SI.list[[2]]$values
subclusters.ast$Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(subclusters.ast$Subclusters, hypo.integrated@meta.data$integrated_Subcluster)]

subclusters.ast$Cluster <- factor(subclusters.ast$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))


###### Can I calculate how many paralog differences (between zeb and mexican tetra) are shared between cave and surface?
###### If I provide a number, say 80% or something, then I can conclude that the majority of them are shared
###### Just find the % of ast that are also in con.ast

calcParalog2 <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, species.3 = species.3, mart.1 = mart[[1]], mart.2 = mart[[2]], ngenes.1 = 32191, ngenes.2 = 25271, i = i) {
  names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2), names(species.3)))
  conserved <- conserved[names]
  species.1 <- species.1[names]
  species.2 <- species.2[names]
  species.3 <- species.3[names]
  
  genes.1 <- setdiff(row.names(species.1[[i]]), row.names(conserved[[i]]))
  genes.2 <- setdiff(row.names(species.2[[i]]), row.names(conserved[[i]]))
  genes.3 <- setdiff(row.names(species.3[[i]]), row.names(conserved[[i]]))
  
  paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(row.names(conserved[[i]]), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(row.names(conserved[[i]]), mart.2$Gene.name)])
  paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart.1$Gene.name)])
  paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart.2$Gene.name)])
  paralog.3 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(genes.3, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(genes.3, mart.2$Gene.name)])
  paralog.union <- union(paralog.con, union(paralog.1, paralog.2))
  
  paralog.1 <- paralog.1[paralog.1 %in% mart.2$Gene.name]
  paralog.2 <- paralog.2[paralog.2 %in% mart.1$Gene.name]
  paralog.3 <- paralog.3[paralog.3 %in% mart.1$Gene.name]
  
  a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.2)])
  b1 <- length(genes.1) - a1
  c1 <- length(union(paralog.con, paralog.2)) - a1
  d1 <- ngenes.1 - b1
  
  a2 <- length(genes.2[genes.2 %in% union(paralog.con, paralog.1)])
  b2 <- length(genes.2) - a2
  c2 <- length(union(paralog.con, paralog.1)) - a2
  d2 <- ngenes.2 - b2
  
  vec1 <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2)
  
  a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.3)])
  b1 <- length(genes.1) - a1
  c1 <- length(union(paralog.con, paralog.3)) - a1
  d1 <- ngenes.1 - b1
  
  a2 <- length(genes.3[genes.3 %in% union(paralog.con, paralog.1)])
  b2 <- length(genes.3) - a2
  c2 <- length(union(paralog.con, paralog.1)) - a2
  d2 <- ngenes.2 - b2
  
  a3 <- length(row.names(conserved[[i]]))
  b3 <- length(genes.1)
  c3 <- length(genes.2)
  d3 <- ngenes.1 - sum(b3, c3)
  
  vec2 <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2, a3=a3, b3=b3, c3=c3, d3=d3)
  vec <- c(vec1, vec2)
  return(vec)
}

paralog.numbers <- lapply(seq_along(gene.lists.pos[[1]]), function(x) calcParalog2(conserved = gene.lists.pos[["cluster.conserved"]], species.1 = gene.lists.pos[["cluster.zebrafish"]], species.2 = gene.lists.pos[["cluster.astyanax"]], species.3 = gene.lists.pos[["cluster.conserved.ast"]], i = x))


clusters.comp <- data.frame(do.call(rbind, paralog.numbers))
row.names(clusters.comp) <- names(gene.lists.pos[[1]])
clusters.comp$percent.para.1 <- unlist(apply(clusters.comp, 1, function(x) x[1]/sum(x[1], x[3])))*100
clusters.comp$percent.para.2 <- unlist(apply(clusters.comp, 1, function(x) x[5]/sum(x[5], x[7])))*100
clusters.comp$Clusters <- row.names(clusters.comp)
clusters.comp$similarityindex <- SI.list[[1]]$values

clusters.comp$Clusters <- factor(clusters.comp$Clusters, levels = levels(hypo.integrated@meta.data$integrated_Cluster))


## For subclusters
intersect <- intersect(names(gene.lists.pos$subcluster.zebrafish), names(gene.lists.pos$subcluster.conserved.ast))

paralog.numbers.sub <- lapply(intersect, function(x) calcParalog2(conserved = gene.lists.pos[["subcluster.conserved"]], species.1 = gene.lists.pos[["subcluster.zebrafish"]], species.2 = gene.lists.pos[["subcluster.astyanax"]], species.3 = gene.lists.pos[["subcluster.conserved.ast"]], i = x))


subclusters.comp <- data.frame(do.call(rbind, paralog.numbers.sub))
row.names(subclusters.comp) <- intersect
subclusters.comp$percent.para.1 <- unlist(apply(subclusters.comp, 1, function(x) x[1]/sum(x[1], x[3])))*100
subclusters.comp$percent.para.2 <- unlist(apply(subclusters.comp, 1, function(x) x[5]/sum(x[5], x[7])))*100
subclusters.comp$Subclusters <- row.names(subclusters.comp)
subclusters.comp$similarityindex <- SI.list[[2]]$values[match(subclusters.comp$Subclusters, SI.list[[2]]$Subcluster)]
subclusters.comp$Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(subclusters.comp$Subclusters, hypo.integrated@meta.data$integrated_Subcluster)]

subclusters.comp$Cluster <- factor(subclusters.comp$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))



## Save results

para.results <- list(clusters, subclusters, clusters.ast, subclusters.ast, clusters.comp, subclusters.comp)
names(para.results) <- c("clusters", "subclusters", "clusters.ast", "subclusters.ast", "clusters.comp", "subclusters.comp")

saveRDS(para.results, file = paste("Paralog-results_markers-trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))


