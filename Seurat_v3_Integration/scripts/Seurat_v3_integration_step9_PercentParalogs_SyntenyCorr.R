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
a = 1.5
b = 2
f = 0.1

gene.lists.pos <- readRDS(paste("drift_gene_lists_pos_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

SI.list <- readRDS(paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

## Filter/update SI.list using synteny corrected orthology (in one direction) - correct astyanax

mart.zeb.89 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/ens89_mart_export_GRCz10_orthologs.txt")

orthologs.3 <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/synteny_orthologs_filt.rds")

# Correct the astyanax lists, then re-calculate the conserved lists using the correct astyanax list and zeb lists

idents <- names(gene.lists.pos$cluster.conserved)

for(i in 1:length(idents)) {
  gene.lists.pos$cluster.astyanax[[idents[[i]]]]$row.names <- row.names(gene.lists.pos$cluster.astyanax[[idents[[i]]]])
  gene.lists.pos$cluster.astyanax[[idents[[i]]]]$corrected <- orthologs.3$Zebrafish.gene.name[match(gene.lists.pos$cluster.astyanax[[idents[[i]]]]$row.names, orthologs.3$Cave.fish.gene.name)]
  gene.lists.pos$cluster.astyanax[[idents[[i]]]]$corrected[is.na(gene.lists.pos$cluster.astyanax[[idents[[i]]]]$corrected)] <- row.names(gene.lists.pos$cluster.astyanax[[idents[[i]]]])[is.na(gene.lists.pos$cluster.astyanax[[idents[[i]]]]$corrected)]
  row.names(gene.lists.pos$cluster.astyanax[[idents[[i]]]]) <- make.names(gene.lists.pos$cluster.astyanax[[idents[[i]]]]$corrected, unique = T)
}

idents2 <- names(gene.lists.pos$subcluster.conserved)

for(i in 1:length(idents)) {
  gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$row.names <- row.names(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]])
  gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$corrected <- orthologs.3$Zebrafish.gene.name[match(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$row.names, orthologs.3$Cave.fish.gene.name)]
  gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$corrected[is.na(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$corrected)] <- row.names(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]])[is.na(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$corrected)]
  row.names(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]) <- make.names(gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]$corrected, unique = T)
}

for(i in 1:length(gene.lists.pos$cluster.conserved)){
  cell_type1 <- gene.lists.pos$cluster.zebrafish[[idents[[i]]]]
  cell_type2 <- gene.lists.pos$cluster.astyanax[[idents[[i]]]]
  colnames(cell_type1) <- paste("cell_type1", colnames(cell_type1), sep = "_")
  colnames(cell_type2) <- paste("cell_type2", colnames(cell_type2), sep = "_")
  index <- intersect(row.names(cell_type1), row.names(cell_type2))
  df <- cbind(cell_type1[index,], cell_type2[index,])
  df$minimup_p_val <- apply(df[,grep("p_val$", colnames(df))], 1, function(x) metap::minimump(x)$p)
  df <- df[df$minimup_p_val < 0.05,]
  gene.lists.pos$cluster.conserved[[i]] <- df
}
names(gene.lists.pos$cluster.conserved) <- idents

for(i in 1:length(gene.lists.pos$subcluster.conserved)){
  cell_type1 <- gene.lists.pos$subcluster.zebrafish[[idents2[[i]]]]
  cell_type2 <- gene.lists.pos$subcluster.astyanax[[idents2[[i]]]]
  colnames(cell_type1) <- paste("cell_type1", colnames(cell_type1), sep = "_")
  colnames(cell_type2) <- paste("cell_type2", colnames(cell_type2), sep = "_")
  index <- intersect(row.names(cell_type1), row.names(cell_type2))
  df <- cbind(cell_type1[index,], cell_type2[index,])
  df$minimup_p_val <- apply(df[,grep("p_val$", colnames(df))], 1, function(x) metap::minimump(x)$p)
  df <- df[df$minimup_p_val < 0.05,]
  gene.lists.pos$subcluster.conserved[[i]] <- df
}
names(gene.lists.pos$subcluster.conserved) <- idents2



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
	con.genes <- row.names(conserved[[i]])
	
	paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(con.genes, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(con.genes, mart.2$Gene.name)])
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

subtypes <- data.frame(do.call(rbind, paralog.numbers))
row.names(subtypes) <- names(gene.lists.pos[[1]])
subtypes$percent.para.1 <- unlist(apply(subtypes, 1, function(x) x[1]/sum(x[1], x[3])))*100
subtypes$percent.para.2 <- unlist(apply(subtypes, 1, function(x) x[5]/sum(x[5], x[7])))*100
# subtypes$fisher.pval <- as.numeric(fisher.results[,1])
subtypes$fisher.odds.1 <- as.numeric(fisher.results.1[,3])
subtypes$fisher.pvalue.1 <- as.numeric(fisher.results.1[,1])
subtypes$fisher.odds.2 <- as.numeric(fisher.results.2[,3])
subtypes$fisher.pvalue.2 <- as.numeric(fisher.results.2[,1])
subtypes$Clusters <- row.names(subtypes)
subtypes$similarityindex <- SI.list[[1]]$values

subtypes$Clusters <- factor(subtypes$Clusters, levels = levels(hypo.integrated@meta.data$integrated_Cluster))


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
subclusters$similarityindex <- SI.list[[2]]$values
subclusters$Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(subclusters$Subclusters, hypo.integrated@meta.data$integrated_Subcluster)]

subclusters$Cluster <- factor(subclusters$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

## Save results

para.results <- list(subtypes, subclusters)

saveRDS(para.results, file = "Paralog_results.rds")

