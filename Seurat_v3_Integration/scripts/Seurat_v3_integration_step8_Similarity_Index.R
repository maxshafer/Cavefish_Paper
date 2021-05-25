library(dplyr)
library(tidyr)
library(Seurat)
library(tictoc)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")
Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

## Load trinarized filtered marker gene lists

a = 1.5
b = 2
f = 0.1

gene.lists.pos <- readRDS(paste("drift_gene_lists_pos_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

## Calculate the drift index!
## Takes 3 lists of marker genes, conserved, speices.1 and species.2
## Takes the intersection of all 3 Cluster names, and subset/order lists by that
## Then run Drift Index caculation for each Cluster
## Maybe add a part for using a cutoff (p-val, logfc, percent)
## Add a part for putting number of genes

calcDriftIndex <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	SI <- list()
	if (is.null(subset)) {
		SI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (nrow(conserved[[x]]) / nrow(species.1[[x]]))) * (1 - (nrow(conserved[[x]]) / nrow(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		SI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		SI <- lapply(seq_along(conserved), function(x) 1 - sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
	}
	}
	names(SI) <- names
	return(SI)
}

## For Cluster
drift.index.Cluster <- calcDriftIndex(conserved = gene.lists.pos$cluster.conserved, species.1 = gene.lists.pos$cluster.zebrafish, species.2 = gene.lists.pos$cluster.astyanax)
SI <- data.frame(Cluster = names(drift.index.Cluster), values = unlist(drift.index.Cluster))

drift.index.Cluster.ast <- calcDriftIndex(conserved = gene.lists.pos$cluster.conserved.ast, species.1 = gene.lists.pos$cluster.surface, species.2 = gene.lists.pos$cluster.cave)
SI.ast <- data.frame(Cluster = c(names(unlist(drift.index.Cluster.ast))), values = c(unlist(drift.index.Cluster.ast)))

## For Subclusters
drift.index.Subcluster <- calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax)
SI.sub <- data.frame(Subcluster = c(names(unlist(drift.index.Subcluster))), values = c(unlist(drift.index.Subcluster)))

drift.index.Subcluster.ast <- calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved.ast, species.1 = gene.lists.pos$subcluster.surface, species.2 = gene.lists.pos$subcluster.cave)
SI.ast.sub <- data.frame(Subcluster = c(names(unlist(drift.index.Subcluster.ast))), values = c(unlist(drift.index.Subcluster.ast)))

# OK, add Cluster as a column to subcluster df

index <- unique(data.frame(Cluster = hypo.integrated@meta.data$integrated_Cluster, Subcluster = hypo.integrated@meta.data$integrated_Subcluster))
SI.sub$Cluster <- index$Cluster[match(SI.sub$Subcluster, index$Subcluster)]
SI.sub$Cluster <- factor(SI.sub$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

index.ast <- unique(data.frame(Cluster = hypo.integrated.ast@meta.data$integrated_Cluster, Subcluster = hypo.integrated.ast@meta.data$integrated_Subcluster))
SI.ast.sub$Cluster <- index.ast$Cluster[match(SI.ast.sub$Subcluster, index.ast$Subcluster)]
SI.ast.sub$Cluster <- factor(SI.ast.sub$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

# Caculate drift for only TFs, NPS/NTS etc
# Get Go lists and use only those that are marker genes

go_lists <- list.files("../Seurat_v3_Integration/SCENIC")[grep("GO", list.files("../Seurat_v3_Integration/SCENIC"))]
go_lists <- lapply(paste("../Seurat_v3_Integration/SCENIC/", go_lists, sep = ""), function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files("../Seurat_v3_Integration/SCENIC")[grep("GO", list.files("../Seurat_v3_Integration/SCENIC"))]
names(go_lists)
go_lists[[14]] <- Reduce(union, list(go_lists[["GO_dopamine.csv"]], go_lists[["GO_epinephrine.csv"]], go_lists[["GO_GABA.csv"]], go_lists[["GO_glutamine.csv"]], go_lists[["GO_neuropeptide.csv"]]))

## Caculate SI for various sets of genes (TFs, NPs, NTs)
## Should probably do this as a 'corrected' SI, will show that the paralogs don't contribute to TF divergence
## Then decide if 'all' should be all genes, or non-TF genes (invert)

SI2 <- data.frame(cell_type = c(names(gene.lists.pos$subcluster.conserved)), 
                  all = c(unlist(calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax, invert = F))), 
                  TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax, subset = go_lists[["GO_TF_list.csv"]], invert = F))), 
                  NP_NTS = c(unlist(calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax, subset = go_lists[[14]], invert = F))))

# SI2 <- SI2[!is.na(SI2$NP_NTS),]
SI.sub.GO <- reshape2::melt(SI2)
SI.sub.GO$variable <- factor(SI.sub.GO$variable, levels = c("all", "NP_NTS", "TFs"))
SI.sub.GO$value[is.infinite(SI.sub.GO$value)] <- 0
# Remove rows with NA values (for non-neuronal NPS_NTS)
SI.sub.GO <- SI.sub.GO[!is.na(SI.sub.GO$value),]

## Save

SI.list <- list(SI, SI.sub, SI.ast, SI.ast.sub, SI.sub.GO)
names(SI.list) <- c("SI", "SI.sub", "SI.ast", "SI.ast.sub", "SI.sub.GO")

# saveRDS(SI.list, file = paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

## Calculate Drift between all cell types both within and between species

## OK, conserved marker genes are found by comparing the p vals between lists
## can use metap::minimump to exclude genes which appear in both gene lists, but are not significant
## compares the two p_val (not adj_p_val), and outputs a minimump_p_val
## OR just the intersect, or calculate the intersect, and then also the minimump
## so need to make data.frames for all cell types comparisons? and then compute the minimump

## Populate a matrix with the SI number, using the row and col names

calcDriftIndex2 <- function(cell_type1 = cell_type1, cell_type2 = cell_type2) {
  colnames(cell_type1) <- paste("cell_type1", colnames(cell_type1), sep = "_")
  colnames(cell_type2) <- paste("cell_type2", colnames(cell_type2), sep = "_")
  index <- intersect(row.names(cell_type1), row.names(cell_type2))
  df <- cbind(cell_type1[index,], cell_type2[index,])
  df$minimup_p_val <- apply(df[,grep("p_val$", colnames(df))], 1, function(x) metap::minimump(x)$p)
  df <- df[df$minimup_p_val < 0.05,]
  
  SI <- 1 - sqrt( abs( (1 - (nrow(df) / nrow(cell_type1))) * (1 - (nrow(df) / nrow(cell_type2))) ) )
  return(SI)
}

matrix.zeb <- list()
for (i in 1:length(gene.lists.pos$cluster.zebrafish)) {
  matrix.zeb[[i]] <- lapply(gene.lists.pos$cluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$cluster.zebrafish[[i]], cell_type2 = x))
}
SI.list[[6]] <- matrix.zeb

matrix.ast <- list()
for (i in 1:length(gene.lists.pos$cluster.astyanax)) {
  matrix.ast[[i]] <- lapply(gene.lists.pos$cluster.astyanax, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$cluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[7]] <- matrix.ast

matrix.zeb.sub <- list()
for (i in 1:length(gene.lists.pos$subcluster.zebrafish)) {
  matrix.zeb.sub[[i]] <- lapply(gene.lists.pos$subcluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$subcluster.zebrafish[[i]], cell_type2 = x))
}
SI.list[[8]] <- matrix.zeb.sub

matrix.ast.sub <- list()
for (i in 1:length(gene.lists.pos$subcluster.astyanax)) {
  matrix.ast.sub[[i]] <- lapply(gene.lists.pos$subcluster.astyanax, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$subcluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[9]] <- matrix.ast.sub

matrix.ast2zeb <- list()
for (i in 1:length(gene.lists.pos$cluster.astyanax)) {
  matrix.ast2zeb[[i]] <- lapply(gene.lists.pos$cluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$cluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[10]] <- matrix.ast2zeb

matrix.ast2zeb.sub <- list()
for (i in 1:length(gene.lists.pos$subcluster.astyanax)) {
  matrix.ast2zeb.sub[[i]] <- lapply(gene.lists.pos$subcluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = gene.lists.pos$subcluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[11]] <- matrix.ast2zeb.sub

names(SI.list) <- c("SI", "SI.sub", "SI.ast", "SI.ast.sub", "SI.sub.GO", "SI.zeb2zeb", "SI.ast2ast", "SI.zeb2zeb.sub", "SI.ast2ast.sub", "SI.ast2zeb", "SI.ast2zeb.sub")

# saveRDS(SI.list, file = paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))


## Make corrected SI
# Subtract paralog numbers from species.1 and species.2 values, and add them to conserved, then calculate SI

SI.list <- readRDS(paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))

mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")


calcCorrDriftIndex <- function(cell_type1 = cell_type1, cell_type2 = cell_type2, mart.1 = mart[[1]], mart.2 = mart[[2]], subset = NULL) {
  colnames(cell_type1) <- paste("cell_type1", colnames(cell_type1), sep = "_")
  colnames(cell_type2) <- paste("cell_type2", colnames(cell_type2), sep = "_")
  
  if (!(is.null(subset))) {
    cell_type1 <- cell_type1[row.names(cell_type1) %in% subset,]
    cell_type2 <- cell_type2[row.names(cell_type2) %in% subset,]
  }
  
  index <- intersect(row.names(cell_type1), row.names(cell_type2))
  df <- cbind(cell_type1[index,], cell_type2[index,])
  df$minimup_p_val <- apply(df[,grep("p_val$", colnames(df))], 1, function(x) metap::minimump(x)$p)
  df <- df[df$minimup_p_val < 0.05,]
  
  paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[chmatch(row.names(df), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[chmatch(row.names(df), mart.2$Gene.name)])
  paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[chmatch(row.names(cell_type2), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[chmatch(row.names(cell_type2), mart.2$Gene.name)])
  paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[chmatch(row.names(cell_type1), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[chmatch(row.names(cell_type1), mart.2$Gene.name)])
  
  paralog.1 <- paralog.1[paralog.1 %chin% mart.2$Gene.name]
  paralog.2 <- paralog.2[paralog.2 %chin% mart.1$Gene.name]
  
  n.genes.1 <- cell_type1[!(cell_type1 %chin% index)] # species-specific genes
  n.genes.2 <- cell_type2[!(cell_type2 %chin% index)]
  
  a1 <- length(row.names(n.genes.1)[row.names(n.genes.1) %chin% union(paralog.con, paralog.2)]) # n.genes.1 that are paralogs of a conserved or species 2 gene (should be subtracted from cell_type1 and added to con)
  a2 <- length(row.names(n.genes.2)[row.names(n.genes.2) %chin% union(paralog.con, paralog.1)])
  
  n.genes.1 <- nrow(n.genes.1) - a1
  n.genes.2 <- nrow(n.genes.2) - a2
  
  n.con <- length(row.names(df)) + a1 + a2
  
  SI.corr <- 1- sqrt( abs( (1 - (n.con / n.genes.1)) * (1 - (n.con / n.genes.2)) ) )
  SI <- 1 - sqrt( abs( (1 - (nrow(df) / nrow(cell_type1))) * (1 - (nrow(df) / nrow(cell_type2))) ) )
  return(SI.corr)
}


matrix.ast2zeb.sub.corr <- list()
for (i in 1:length(gene.lists.pos$subcluster.astyanax)) {
  tic(paste("finished", names(trinarized.exp$subcluster.astyanax)[[i]], "in", sep = " "))
  matrix.ast2zeb.sub.corr[[i]] <- lapply(gene.lists.pos$subcluster.zebrafish, function(x) calcCorrDriftIndex(cell_type1 = gene.lists.pos$subcluster.astyanax[[i]], cell_type2 = x))
  toc()
}
names(matrix.ast2zeb.sub.corr) <- names(gene.lists.pos$subcluster.astyanax)

SI.list[[12]] <- matrix.ast2zeb.sub.corr

## Save

names(SI.list) <- c("SI", "SI.sub", "SI.ast", "SI.ast.sub", "SI.sub.GO", "SI.zeb2zeb", "SI.ast2ast", "SI.zeb2zeb.sub", "SI.ast2ast.sub", "SI.ast2zeb", "SI.ast2zeb.sub", "SI.ast2zeb.sub.corr")

saveRDS(SI.list, file = paste("SI_results_trinarized_a", a, "_b", b, "_f",f,".rds", sep = ""))



## Calculate corrected SI, based on gene sets
## cell_type1 and cell_type2 are the dfs for each, so need to lapply/for looop here as well
all <- lapply(seq_along(gene.lists.pos$subcluster.zebrafish), function(x) calcCorrDriftIndex(cell_type1 = gene.lists.pos$subcluster.zebrafish[[x]], cell_type2 = gene.lists.pos$subcluster.astyanax[[x]]))

TFs <- lapply(seq_along(gene.lists.pos$subcluster.zebrafish), function(x) calcCorrDriftIndex(cell_type1 = gene.lists.pos$subcluster.zebrafish[[x]], cell_type2 = gene.lists.pos$subcluster.astyanax[[x]], subset = go_lists[["GO_TF_list.csv"]]))

NP_NTS <- lapply(seq_along(gene.lists.pos$subcluster.zebrafish), function(x) calcCorrDriftIndex(cell_type1 = gene.lists.pos$subcluster.zebrafish[[x]], cell_type2 = gene.lists.pos$subcluster.astyanax[[x]], subset = go_lists[[14]]))

SI2 <- data.frame(cell_type = c(names(gene.lists.pos$subcluster.zebrafish)), all = unlist(all), TFs = unlist(TFs), NP_NTS = unlist(NP_NTS))

SI2 <- data.frame(cell_type = c(names(gene.lists.pos$subcluster.conserved)), 
                  all = c(unlist(calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax, invert = F))), 
                  TFs = c(unlist(calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax, subset = go_lists[["GO_TF_list.csv"]], invert = F))), 
                  NP_NTS = c(unlist(calcDriftIndex(conserved = gene.lists.pos$subcluster.conserved, species.1 = gene.lists.pos$subcluster.zebrafish, species.2 = gene.lists.pos$subcluster.astyanax, subset = go_lists[[14]], invert = F))))

# SI2 <- SI2[!is.na(SI2$NP_NTS),]
SI.sub.GO.corr <- reshape2::melt(SI2)
SI.sub.GO.corr$variable <- factor(SI.sub.GO.corr$variable, levels = c("all", "NP_NTS", "TFs"))
SI.sub.GO.corr$value[is.infinite(SI.sub.GO.corr$value)] <- 0
# Remove rows with NA values (for non-neuronal NPS_NTS)
SI.sub.GO.corr <- SI.sub.GO.corr[!is.na(SI.sub.GO.corr$value),]




# ## Try with OG clustering
# 
# matrix.ast2zeb.og <- list()
# for (i in 1:length(ast.markers[[1]])) {
#   matrix.ast2zeb.og[[i]] <- lapply(zeb.markers[[1]], function(x) calcDriftIndex2(cell_type1 = ast.markers[[1]][[i]], cell_type2 = x))
# }
# 
# names(matrix.ast2zeb.og) <- names(ast.markers[[1]])
# matrix.ast2zeb.og.2 <- reshape2::melt(matrix.ast2zeb.og)
# ggplot(matrix.ast2zeb.og.2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_viridis(direction = -1)
# 
# 
# names(matrix.ast2zeb.sub) <- names(ast.markers[[1]])
# test <- lapply(matrix.ast2zeb.sub, function(x) unlist(x))
# test <- lapply(test, function(x) 1-x)
# test2 <- lapply(test, function(x) x/max(x))
# test2 <- reshape2::melt(test2)
# test2$L2 <- rep(names(matrix.ast2zeb.sub[[1]]),151)
# test2$L1 <- factor(test2$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# test2$L2 <- factor(test2$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# ggplot(test2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8)) + scale_fill_viridis(direction = 1)
# 
# test2 <- matrix.ast2zeb.og.2 %>% group_by(L1) %>% top_n(-1, value)
