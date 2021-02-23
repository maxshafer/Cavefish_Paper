library(dplyr)
library(fdrtool)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

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

gene.lists.pos <- readRDS(file = "drift_gene_lists_pos_2.rds")
SI.list <- readRDS(file = "SI_results_trinarized_markers.rds")

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


paralog.numbers <- lapply(seq_along(gene.lists.pos[[1]]), function(x) calcParalog(conserved = gene.lists.pos[[1]], species.1 = gene.lists.pos[[2]], species.2 = gene.lists.pos[[3]], i = x))

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
clusters$driftindex <- SI.list[[1]]$values

clusters$Clusters <- factor(clusters$Clusters, levels = levels(hypo.integrated@meta.data$integrated_Cluster))


#### For subclusters!

paralog.numbers.sub <- lapply(seq_along(gene.lists.pos[[7]]), function(x) calcParalog(conserved = gene.lists.pos[[7]], species.1 = gene.lists.pos[[8]], species.2 = gene.lists.pos[[9]], i = x))

fisher.results.1 <- lapply(paralog.numbers.sub, function(x) fisher.test(matrix(data = x[1:4], ncol = 2, nrow = 2)))
fisher.results.1 <- do.call(rbind, fisher.results.1)
row.names(fisher.results.1) <- names(gene.lists.pos[[7]])
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
subclusters$driftindex <- SI.list[[2]]$values
subclusters$Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(subclusters$Subclusters, hypo.integrated@meta.data$integrated_Subcluster)]

subclusters$Cluster <- factor(subclusters$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

## Save results

para.results <- list(clusters, subclusters)

saveRDS(para.results, file = "Paralog_results_trinarized_markers.rds")

