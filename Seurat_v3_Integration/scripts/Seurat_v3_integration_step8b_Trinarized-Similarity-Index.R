library(dplyr)
library(tidyr)
library(Seurat)
library(tictoc)

# # Set hyperparameters from arguments or in file
# args <- commandArgs(trailingOnly = TRUE)
# 
# print(args)
# 
# a <- as.numeric(args[[1]])
# b <- as.numeric(args[[2]])
# f <- as.numeric(args[[3]])

a = 1.5
b = 2
f = 0.1

# Load files

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")
Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

# Load trinarized gene lists

trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# # Load hk filtered file
# trinarized.exp <- readRDS(file = paste("trinarized_expression_hk-filtered_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

## Calculate the similarity index!
## Takes 2 lists of marker genes, speices.1 and species.2
## Takes the intersection of both Cluster names, and subset/order lists by that
## Then run Similarity Index caculation for each Cluster

calcDriftIndex <- function(species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
  names <- Reduce(intersect, list(names(species.1), names(species.2)))
  
  conserved <- lapply(names, function(x) intersect(names(species.1[[x]]), names(species.2[[x]])))
  species.1 <- lapply(names, function(x) names(species.1[[x]]))
  species.2 <- lapply(names, function(x) names(species.2[[x]]))
  
  DI <- list()
  if (is.null(subset)) {
    DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( (1 - (length(conserved[[x]]) / length(species.1[[x]]))) * (1 - (length(conserved[[x]]) / length(species.2[[x]]))) ))
  } else {
    if (invert == FALSE) {
      DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( (1 - (length(conserved[[x]][conserved[[x]] %in% subset]) / length(species.1[[x]][species.1[[x]] %in% subset]))) * (1 - (length(conserved[[x]][conserved[[x]] %in% subset]) / length(species.2[[x]][species.2[[x]] %in% subset])) )))
    } else {
      DI <- lapply(seq_along(conserved), function(x) 1 - sqrt( (1 - (length(conserved[[x]][!(conserved[[x]] %in% subset)]) / length(species.1[[x]][!(species.1[[x]] %in% subset)]))) * (1 - (length(conserved[[x]][!(conserved[[x]] %in% subset)]) / length(species.2[[x]][!(species.2[[x]] %in% subset)])) ) ))
    }
  }
  names(DI) <- names
  return(DI)
}

## For Clusters
trin.SI.Clusters <- calcDriftIndex(species.1 = trinarized.exp$cluster.zebrafish, species.2 = trinarized.exp$cluster.astyanax)
SI <- data.frame(Cluster = names(trin.SI.Clusters), values = unlist(trin.SI.Clusters))

trin.SI.Clusters.ast <- calcDriftIndex(species.1 = trinarized.exp$cluster.zebrafish, species.2 = trinarized.exp$cluster.astyanax)
SI.ast <- data.frame(Cluster = c(names(unlist(trin.SI.Clusters.ast))), values = c(unlist(trin.SI.Clusters.ast)))

## For Subclusters
trin.SI.Subclusters <- calcDriftIndex(species.1 = trinarized.exp$subcluster.zebrafish, species.2 = trinarized.exp$subcluster.astyanax)
SI.sub <- data.frame(Subcluster = c(names(unlist(trin.SI.Subclusters))), values = c(unlist(trin.SI.Subclusters)))

trin.SI.Subclusters.ast <- calcDriftIndex(species.1 = trinarized.exp$subcluster.zebrafish, species.2 = trinarized.exp$subcluster.astyanax)
SI.ast.sub <- data.frame(Subcluster = c(names(unlist(trin.SI.Subclusters.ast))), values = c(unlist(trin.SI.Subclusters.ast)))

# OK, add Cluster as a column to subcluster df

index <- unique(data.frame(Cluster = hypo.integrated@meta.data$integrated_Cluster, Subcluster = hypo.integrated@meta.data$integrated_Subcluster))
SI.sub$Cluster <- index$Cluster[match(SI.sub$Subcluster, index$Subcluster)]
SI.sub$Cluster <- factor(SI.sub$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))

index.ast <- unique(tibble(Cluster = hypo.integrated.ast@meta.data$integrated_Cluster, Subcluster = hypo.integrated.ast@meta.data$integrated_Subcluster))
SI.ast.sub$Cluster <- index.ast$Cluster[match(SI.ast.sub$Subcluster, index.ast$Subcluster)]
SI.ast.sub$Cluster <- factor(SI.ast.sub$Cluster, levels = levels(hypo.integrated@meta.data$integrated_Cluster))



# Caculate similarity for only TFs, NPS/NTS etc
# Get Go lists and use only those that are marker genes

go_lists <- list.files("../Seurat_v3_Integration/SCENIC")[grep("GO", list.files("../Seurat_v3_Integration/SCENIC"))]
go_lists <- lapply(paste("../Seurat_v3_Integration/SCENIC/", go_lists, sep = ""), function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files("../Seurat_v3_Integration/SCENIC")[grep("GO", list.files("../Seurat_v3_Integration/SCENIC"))]
names(go_lists)
go_lists[[14]] <- Reduce(union, list(go_lists[[5]], go_lists[[6]], go_lists[[9]]))

## Plot Subtypes

SI2 <- tibble(cell_type = c(names(trinarized.exp$subcluster.zebrafish)), 
              all = c(unlist(calcDriftIndex(species.1 = trinarized.exp$subcluster.zebrafish, species.2 = trinarized.exp$subcluster.astyanax, invert = F))), 
              TFs = c(unlist(calcDriftIndex(species.1 = trinarized.exp$subcluster.zebrafish, species.2 = trinarized.exp$subcluster.astyanax, subset = go_lists[["GO_TF_list.csv"]], invert = F))), 
              NP_NTS = c(unlist(calcDriftIndex(species.1 = trinarized.exp$subcluster.zebrafish, species.2 = trinarized.exp$subcluster.astyanax, subset = go_lists[[""]], invert = F))))
SI.sub.GO <- reshape2::melt(SI2)
SI.sub.GO$variable <- factor(SI.sub.GO$variable, levels = c("all", "NP_NTS", "TFs"))
SI.sub.GO <- SI.sub.GO[!is.na(SI.sub.GO$variable),]
# SI.sub.GO$value[is.infinite(SI.sub.GO$value)] <- 0


## Save

SI.list <- list(SI, SI.sub, SI.ast, SI.ast.sub, SI.sub.GO)
names(SI.list) <- c("SI", "SI.sub", "SI.ast", "SI.ast.sub", "SI.sub.GO")

# saveRDS(SI.list, file = paste("SI_trinarized_a",a,"_b",b, "_f",f,".rds", sep = ""))


## Calculate Similarity between all cell types both within and between species

## Very easy, just lapply across the lists
## Populate a matrix with the SI number, using the row and col names

calcDriftIndex2 <- function(cell_type1 = cell_type1, cell_type2 = cell_type2) {
  cell_type1 <- names(cell_type1)
  cell_type2 <- names(cell_type2)
  index <- intersect(cell_type1, cell_type2)
  
  DI <- 1 - sqrt( abs( (1 - (length(index) / length(cell_type1))) * (1 - (length(df) / length(cell_type2))) ) )
  return(DI)
}

matrix.zeb <- list()
for (i in 1:length(trinarized.exp$cluster.zebrafish)) {
  matrix.zeb[[i]] <- lapply(trinarized.exp$cluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = trinarized.exp$cluster.zebrafish[[i]], cell_type2 = x))
}
SI.list[[6]] <- matrix.zeb

matrix.ast <- list()
for (i in 1:length(trinarized.exp$cluster.astyanax)) {
  matrix.ast[[i]] <- lapply(trinarized.exp$cluster.astyanax, function(x) calcDriftIndex2(cell_type1 = trinarized.exp$cluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[7]] <- matrix.ast

matrix.zeb.sub <- list()
for (i in 1:length(trinarized.exp$subcluster.zebrafish)) {
  matrix.zeb.sub[[i]] <- lapply(trinarized.exp$subcluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = trinarized.exp$subcluster.zebrafish[[i]], cell_type2 = x))
}
SI.list[[8]] <- matrix.zeb.sub

matrix.ast.sub <- list()
for (i in 1:length(trinarized.exp$subcluster.astyanax)) {
  matrix.ast.sub[[i]] <- lapply(trinarized.exp$subcluster.astyanax, function(x) calcDriftIndex2(cell_type1 = trinarized.exp$subcluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[9]] <- matrix.ast.sub

matrix.ast2zeb <- list()
for (i in 1:length(trinarized.exp$cluster.astyanax)) {
  matrix.ast2zeb[[i]] <- lapply(trinarized.exp$cluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = trinarized.exp$cluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[10]] <- matrix.ast2zeb

matrix.ast2zeb.sub <- list()
for (i in 1:length(trinarized.exp$subcluster.astyanax)) {
  matrix.ast2zeb.sub[[i]] <- lapply(trinarized.exp$subcluster.zebrafish, function(x) calcDriftIndex2(cell_type1 = trinarized.exp$subcluster.astyanax[[i]], cell_type2 = x))
}
SI.list[[11]] <- matrix.ast2zeb.sub

names(SI.list) <- c("SI", "SI.sub", "SI.ast", "SI.ast.sub", "SI.sub.GO", "SI.zeb2zeb", "SI.ast2ast", "SI.zeb2zeb.sub", "SI.ast2ast.sub", "SI.ast2zeb", "SI.ast2zeb.sub")

saveRDS(SI.list, file = paste("SI_trinarized_a",a,"_b",b, "_f",f,".rds", sep = ""))
# saveRDS(SI.list, file = paste("SI_trinarized_hk-filtered_a",a,"_b",b, "_f",f,".rds", sep = ""))

## Make corrected SI
# Subtract paralog numbers from species.1 and species.2 values, and add them to conserved, then calculate SI

mart <- list()
mart[[1]] <- read.csv("mart_export_GRCz10_paralogs.txt", head = TRUE)
mart[[2]] <- read.csv("mart_export_AstMex102_paralogs.txt", head = TRUE)
names(mart) <- c("zebrafish", "astyanax")


calcCorrDriftIndex <- function(cell_type1 = cell_type1, cell_type2 = cell_type2, mart.1 = mart[[1]], mart.2 = mart[[2]]) {
  cell_type1 <- names(cell_type1)
  cell_type2 <- names(cell_type2)
  
  index <- intersect(cell_type1, cell_type2)
  
  paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[chmatch(index, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[chmatch(index, mart.2$Gene.name)])
  paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[chmatch(cell_type2, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[chmatch(cell_type2, mart.2$Gene.name)])
  paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[chmatch(cell_type1, mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[chmatch(cell_type1, mart.2$Gene.name)])
  
  # paralog.1 <- paralog.1[paralog.1 %in% mart.2$Gene.name]
  # paralog.2 <- paralog.2[paralog.2 %in% mart.1$Gene.name]
  
  paralog.1 <- paralog.1[paralog.1 %chin% mart.2$Gene.name]
  paralog.2 <- paralog.2[paralog.2 %chin%  mart.1$Gene.name]
  
  n.genes.1 <- cell_type1[!(cell_type1 %chin% index)] # species-specific genes
  n.genes.2 <- cell_type2[!(cell_type2 %chin% index)]
  
  a1 <- length(n.genes.1[n.genes.1 %chin% union(paralog.con, paralog.2)]) # n.genes.1 that are paralogs of a conserved or species 2 gene (should be subtracted from cell_type1 and added to con)
  a2 <- length(n.genes.2[n.genes.2 %chin% union(paralog.con, paralog.1)])
  
  n.genes.1 <- length(cell_type1) - a1
  n.genes.2 <- length(cell_type2) - a2
  
  n.con <- length(index) + a1 + a2
  
  DI.corr <- 1- sqrt( abs( (1 - (n.con / n.genes.1)) * (1 - (n.con / n.genes.2)) ) )
  return(DI.corr)
}


matrix.ast2zeb.sub.corr <- list()
for (i in 1:length(trinarized.exp$subcluster.astyanax)) {
  tic(paste("finished", names(trinarized.exp$subcluster.astyanax)[[i]], "in", sep = " "))
  matrix.ast2zeb.sub.corr[[i]] <- lapply(trinarized.exp$subcluster.zebrafish, function(x) calcCorrDriftIndex(cell_type1 = trinarized.exp$subcluster.astyanax[[i]], cell_type2 = x))
  toc()
}
names(matrix.ast2zeb.sub.corr) <- names(trinarized.exp$subcluster.astyanax)

SI.list[[12]] <- matrix.ast2zeb.sub.corr

names(SI.list[[12]]) <- c("SI.ast2zeb.sub.corr")

saveRDS(SI.list, file = paste("SI_trinarized_a",a,"_b",b, "_f",f,".rds", sep = ""))
# saveRDS(SI.list, file = paste("SI_trinarized_hk-filtered_a",a,"_b",b, "_f",f,".rds", sep = ""))


# # a = 1.5
# # b = 2
# # f = 0.35
# # 
# # SI.list <- readRDS(file = paste("SI_trinarized_a",a,"_b",b, "_f",f,".rds", sep = ""))
# # 
# ## Make non row-scaled heatmap
# matrix.ast2zeb.sub <- SI.list[[11]]
# # matrix.ast2zeb.sub <- hk.rm.11.035
# names(matrix.ast2zeb.sub) <- names(matrix.ast2zeb.sub[[1]])
# matrix2 <- reshape2::melt(matrix.ast2zeb.sub)
# matrix2$L1 <- unlist(lapply(names(matrix.ast2zeb.sub), function(x) rep(x, 151)))
# matrix2$L1 <- factor(matrix2$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# matrix2$L2 <- factor(matrix2$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# DI.matrix <- ggplot(matrix2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_fill_viridis(direction = 1)
# 
# ## Make row-scaled heatmap
# test <- lapply(matrix.ast2zeb.sub, function(x) unlist(x))
# # test <- lapply(test, function(x) x[!(names(x) %in% "Neuronal_07_4")]) # removes the zeb subcluster from each of the astyanax subcluster lists
# # test <- test[!(names(test) %in% "Neuronal_07_4")] # this gets rid of the astyanax cluster (row, not column)
# # test <- lapply(test, function(x) 1-x)
# test2 <- lapply(test, function(x) x/max(x))
# test2 <- reshape2::melt(test2)
# test2$L2 <- rep(names(test[[1]]),151)
# # test2 <- test2[!(test2$L2 == "Neuronal_07_4"),] # removes the column, but only after normalisation
# test2$L1 <- factor(test2$L1, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# test2$L2 <- factor(test2$L2, levels = levels(hypo.integrated@meta.data$integrated_Subcluster))
# DI.matrix.scaled <- ggplot(test2, aes(x = L2, y = L1, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 8)) + scale_fill_viridis(direction = 1)
# 
# dev.new()
# DI.matrix + DI.matrix.scaled
# 
# DI.matrix.scaled.035.hk <- DI.matrix.scaled
# # 
# # 
# # plot.list <- DI.plots.01 + DI.plots.02 + DI.plots.035 + plot_layout(nrow = 3)
