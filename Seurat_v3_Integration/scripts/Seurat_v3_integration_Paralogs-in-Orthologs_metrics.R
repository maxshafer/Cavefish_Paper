library(Matrix)
library(dplyr)
library(ggplot2)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(data.table)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

### Load mart files, and determine gene pair lists
# Read in the file from biomart web tool (R package doesn't work, because file is so huge)
mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")

### Make new gene list that is all the subcluster genes, but named by subtype
gene.lists.pos <- readRDS(file = "drift_gene_lists_pos.rds")
str(gene.lists.pos, max.level = 1)

## Subset this list for non-paralogs

# Need to do this for (1) 2,3, and (7) 8,9

gene.lists.para <- gene.lists.pos

subtypes <- names(gene.lists.pos[[1]])
subclusters <- names(gene.lists.pos[[7]])


# Need to make a df with pairs of paralogs which are expressed in orthologous cell types
# Add column with the orthology confidence score for those two genes, plus the orthology score for each gene with it's called ortholog

gene.lists.para[[2]] <- lapply(seq_along(subtypes), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[2]][[subtypes[x]]]), row.names(gene.lists.pos[[1]][[subtypes[x]]])) # species.1 specific genes
  genes.2 <- setdiff(row.names(gene.lists.pos[[3]][[subtypes[x]]]), row.names(gene.lists.pos[[1]][[subtypes[x]]])) # species.2 specific genes
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subtypes[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subtypes[x]]]), mart[[2]]$Gene.name)])
  paralog.2 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart[[2]]$Gene.name)])
  paralog.2 <- union(paralog.con, paralog.2)
  paralog.2 <- paralog.2[paralog.2 %in% mart[[1]]$Gene.name]
  return(gene.lists.pos[[2]][[subtypes[x]]][genes.1[genes.1 %in% paralog.2],])
})
gene.lists.para[[3]] <- lapply(seq_along(subtypes), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[2]][[subtypes[x]]]), row.names(gene.lists.pos[[1]][[subtypes[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[3]][[subtypes[x]]]), row.names(gene.lists.pos[[1]][[subtypes[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subtypes[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subtypes[x]]]), mart[[2]]$Gene.name)])
  paralog.1 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart[[2]]$Gene.name)])
  paralog.1 <- union(paralog.con, paralog.1)
  paralog.1 <- paralog.1[paralog.1 %in% mart[[2]]$Gene.name]
  return(gene.lists.pos[[3]][[subtypes[x]]][genes.2[genes.2 %in% paralog.1],])
})

names(gene.lists.para[[2]]) <- subtypes
names(gene.lists.para[[2]]) <- subtypes

gene.lists.para[[8]] <- lapply(seq_along(subclusters), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[8]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[9]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[2]]$Gene.name)])
  paralog.2 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.2, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.2, mart[[2]]$Gene.name)])
  paralog.2 <- union(paralog.con, paralog.2)
  paralog.2 <- paralog.2[paralog.2 %in% mart[[1]]$Gene.name]
  return(gene.lists.pos[[8]][[subclusters[x]]][genes.1[genes.1 %in% paralog.2],])
})
gene.lists.para[[9]] <- lapply(seq_along(subclusters), function(x) {
  genes.1 <- setdiff(row.names(gene.lists.pos[[8]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  genes.2 <- setdiff(row.names(gene.lists.pos[[9]][[subclusters[x]]]), row.names(gene.lists.pos[[7]][[subclusters[x]]]))
  paralog.con <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(row.names(gene.lists.pos[[1]][[subclusters[x]]]), mart[[2]]$Gene.name)])
  paralog.1 <- union(mart[[1]]$Zebrafish.paralogue.associated.gene.name[match(genes.1, mart[[1]]$Gene.name)], mart[[2]]$Cave.fish.paralogue.associated.gene.name[match(genes.1, mart[[2]]$Gene.name)])
  paralog.1 <- union(paralog.con, paralog.1)
  paralog.1 <- paralog.1[paralog.1 %in% mart[[2]]$Gene.name]
  return(gene.lists.pos[[9]][[subclusters[x]]][genes.2[genes.2 %in% paralog.1],])
})

names(gene.lists.para[[8]]) <- subclusters
names(gene.lists.para[[9]]) <- subclusters

# for loop, or apply over the names of subtypes, and do this for all three gene.lists.pos[[4:6]]
test <- lapply(list(7,8,9), function(y) lapply(names(gene.lists.pos[[1]]), function(x) Reduce(union, lapply(gene.lists.pos[[y]][grep(x, names(gene.lists.pos[[y]]))], function(z) row.names(z)))))
names(test) <- c("conserved.markers.sub.combined", "zebrafish.markers.sub.combined", "astyanax.markers.sub.combined")

for(i in 1:length(test)){
  names(test[[i]]) <- names(gene.lists.pos[[1]])
}

## Append gene lists
gene.lists.pos <- c(gene.lists.pos[1:14], test)


# for loop, or apply over the names of subtypes, and do this for all three gene.lists.pos[[4:6]]
test <- lapply(list(7,8,9), function(y) lapply(names(gene.lists.para[[1]]), function(x) Reduce(union, lapply(gene.lists.para[[y]][grep(x, names(gene.lists.para[[y]]))], function(z) row.names(z)))))
names(test) <- c("conserved.markers.sub.combined", "zebrafish.markers.sub.combined", "astyanax.markers.sub.combined")

for(i in 1:length(test)){
	names(test[[i]]) <- names(gene.lists.pos[[1]])
}

## Append gene lists
gene.lists.para <- c(gene.lists.para[1:14], test)

saveRDS(gene.lists.para, file = "drift_gene_lists_pos_para.rds")


# Load orthology mart files
mart.zeb <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_GRCz10_orthologs.txt")
mart.ast <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_AstMex102_orthologs.txt")

# Compute the mean of Gene order conservation score for paralogs expressed in zebrafish
para <- unlist(lapply(names(gene.lists.para$zebrafish.markers.sub), function(x) {
  gocs <- mart.zeb[mart.zeb$Gene.name %in% row.names(gene.lists.para[["zebrafish.markers.sub"]][[x]]),"Cave.fish.orthology.confidence..0.low..1.high."]
  mean(gocs[!(is.na(gocs))])
}))
markers.con <- unlist(lapply(names(gene.lists.pos$conserved.markers.sub), function(x) {
  gocs <- mart.zeb[mart.zeb$Gene.name %in% row.names(gene.lists.pos[["conserved.markers.sub"]][[x]]),"Cave.fish.orthology.confidence..0.low..1.high."]
  mean(gocs[!(is.na(gocs))])
}))
markers.zeb <- unlist(lapply(names(gene.lists.pos$zebrafish.markers.sub), function(x) {
  gocs <- mart.zeb[mart.zeb$Gene.name %in% row.names(gene.lists.pos[["zebrafish.markers.sub"]][[x]]),"Cave.fish.orthology.confidence..0.low..1.high."]
  mean(gocs[!(is.na(gocs))])
}))

df.conf <- data.frame(markers.con = markers.con, markers.zeb = markers.zeb, para = para)
row.names(df.conf) <- names(gene.lists.para$conserved.markers.sub)
colnames(df.conf) <- c("Conserved marker genes", "Species-specific marker genes", "Species-specific paralogous genes")

orth_conf_metrics <- list()
orth_conf_metrics[[1]] <- df.conf

# Compute the mean of Gene order conservation score for paralogs expressed in zebrafish
para <- unlist(lapply(names(gene.lists.para$zebrafish.markers.sub), function(x) {
  gocs <- mart.zeb[mart.zeb$Gene.name %in% row.names(gene.lists.para[["zebrafish.markers.sub"]][[x]]),"Cave.fish.Gene.order.conservation.score"]
  mean(gocs[!(is.na(gocs))])
}))
markers.con <- unlist(lapply(names(gene.lists.pos$conserved.markers.sub), function(x) {
  gocs <- mart.zeb[mart.zeb$Gene.name %in% row.names(gene.lists.pos[["conserved.markers.sub"]][[x]]),"Cave.fish.Gene.order.conservation.score"]
  mean(gocs[!(is.na(gocs))])
}))
markers.zeb <- unlist(lapply(names(gene.lists.pos$zebrafish.markers.sub), function(x) {
  gocs <- mart.zeb[mart.zeb$Gene.name %in% row.names(gene.lists.pos[["zebrafish.markers.sub"]][[x]]),"Cave.fish.Gene.order.conservation.score"]
  mean(gocs[!(is.na(gocs))])
}))

df.gene.order <- data.frame(markers.con = markers.con, markers.zeb = markers.zeb, para = para)
row.names(df.gene.order) <- names(gene.lists.para$conserved.markers.sub)
colnames(df.gene.order) <- c("Conserved marker genes", "Species-specific marker genes", "Species-specific paralogous genes")

orth_conf_metrics[[2]] <- df.gene.order
saveRDS(orth_conf_metrics, file = "paralogs-orthology-conf-metrics.rds")



