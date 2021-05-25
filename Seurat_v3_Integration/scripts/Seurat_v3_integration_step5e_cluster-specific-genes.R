library(Seurat)
library(Matrix)
library(dplyr)
library(ggplot2)
library(data.table)
library(patchwork)
library(ggrepel)
library(viridis)
library(ggpubr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

a = 1.5
b = 2
f = 0.1

## Load trinarized gene lists

trinarized.genes <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,".rds", sep = ""))

trinarized.genes.2 <- lapply(trinarized.genes, function(x) Reduce(cbind, x))

for(i in 1:length(trinarized.genes.2)) {
  colnames(trinarized.genes.2[[i]]) <- names(trinarized.genes[[i]])
}


trinarized.genes.2 <- lapply(trinarized.genes.2, function(x) lapply(x, function(y) y[y > 0.95]))

trinarized.genes.2 <- lapply(trinarized.genes.2, function(x) ifelse(x > 0.95, 1, 0))

rowsums <- lapply(trinarized.genes.2, function(x) rowSums(x))

trinarized.genes.3 <- lapply(seq_along(trinarized.genes.2), function(x) trinarized.genes.2[[x]][rowsums[[x]] == 1,])
names(trinarized.genes.3) <- names(trinarized.genes.2)

# Calculate the number of cluster specific genes (not expressed in other clusters) for zeb and ast
df.sub <- data.frame(cell_types = colnames(trinarized.genes.3$subcluster.zebrafish), zeb = unlist(lapply(c(1:151), function(x) length(rownames(trinarized.genes.3$subcluster.zebrafish[trinarized.genes.3$subcluster.zebrafish[,x] == 1,])))), ast = unlist(lapply(c(1:151), function(x) length(rownames(trinarized.genes.3$subcluster.astyanax[trinarized.genes.3$subcluster.astyanax[,x] == 1,])))))
df <- data.frame(cell_types = colnames(trinarized.genes.3$cluster.zebrafish), zeb = unlist(lapply(c(1:24), function(x) length(rownames(trinarized.genes.3$cluster.zebrafish[trinarized.genes.3$cluster.zebrafish[,x] == 1,])))), ast = unlist(lapply(c(1:24), function(x) length(rownames(trinarized.genes.3$cluster.astyanax[trinarized.genes.3$cluster.astyanax[,x] == 1,])))))

# Calculate how many of the cluster specific genes are in common between species
df$intersection <- unlist(lapply(c(1:24), function(x) length(intersect(rownames(trinarized.genes.3$cluster.zebrafish[trinarized.genes.3$cluster.zebrafish[,x] == 1,]), rownames(trinarized.genes.3$cluster.astyanax[trinarized.genes.3$cluster.astyanax[,x] == 1,])))))
df.sub$intersection <- unlist(lapply(c(1:151), function(x) length(intersect(rownames(trinarized.genes.3$subcluster.zebrafish[trinarized.genes.3$subcluster.zebrafish[,x] == 1,]), rownames(trinarized.genes.3$subcluster.astyanax[trinarized.genes.3$subcluster.astyanax[,x] == 1,])))))

# Calculate SI?

df$union <- df$zeb + df$ast - 2*df$intersection
df.sub$union <- df.sub$zeb + df.sub$ast - 2*df.sub$intersection

df$SI <- 1 - sqrt( (1 - (df$intersection / df$zeb)) * (1 - (df$intersection / df$ast)) )
df.sub$SI <- 1 - sqrt( (1 - (df.sub$intersection / df.sub$zeb)) * (1 - (df.sub$intersection / df$ast)) )

ggplot(melt(df.sub), aes(x = variable, y = value)) + geom_jitter() + theme_classic() + facet_wrap(~variable, scales = "free", nrow = 1)

DotPlot(hypo.integrated.zeb, features = rownames(trinarized.genes.3$subcluster.zebrafish[trinarized.genes.3$subcluster.zebrafish[,"Neuronal_03_15"] == 1,]), group.by = "integrated_Subcluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

DotPlot(hypo.integrated.zeb, features = rownames(trinarized.genes.3$subcluster.zebrafish[trinarized.genes.3$subcluster.zebrafish[,"Neuronal_03_15"] == 1,]), group.by = "integrated_Subcluster") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))





