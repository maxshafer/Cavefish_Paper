library(dplyr)
library(tidyr)
library(Seurat)

# Load files

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

# Set hyperparameters

a = 1.5
b = 2
f = 0.1

## Load trinarized gene lists

trinarized.genes <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,".rds", sep = ""))

trinarized.exp <- readRDS(file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

# Find housekeeping genes
trinarized.hk <- trinarized.genes
hk <- list()
trinarized.hk$cluster.zebrafish <- Reduce(cbind, trinarized.hk$cluster.zebrafish)
trinarized.hk$cluster.zebrafish <- trinarized.hk$cluster.zebrafish > 0.95
hk[[1]] <- apply(trinarized.hk$cluster.zebrafish, 1, function(x) all(x))[apply(trinarized.hk$cluster.zebrafish, 1, function(x) all(x))]

trinarized.hk$cluster.astyanax <- Reduce(cbind, trinarized.hk$cluster.astyanax)
trinarized.hk$cluster.astyanax <- trinarized.hk$cluster.astyanax > 0.95
hk[[2]] <- apply(trinarized.hk$cluster.astyanax, 1, function(x) all(x))[apply(trinarized.hk$cluster.astyanax, 1, function(x) all(x))]

hk[[1]] <- apply(trinarized.hk$cluster.zebrafish, 1, function(x) all(x))[apply(trinarized.hk$cluster.zebrafish, 1, function(x) length(x[x]) > length(x)*.4)]
hk[[2]] <- apply(trinarized.hk$cluster.astyanax, 1, function(x) all(x))[apply(trinarized.hk$cluster.astyanax, 1, function(x) length(x[x]) > length(x)*.4)]

saveRDS(hk, file = paste("housekeeping_genes_a",a,"_b",b, "_f",f,".rds", sep = ""))

## Calculate housekeeping SI

Gt <- length(intersect(names(hk[[1]][hk[[1]]]), names(hk[[2]][hk[[2]]])))
Ga <- length(names(hk[[1]][hk[[1]]]))
Gb <- length(names(hk[[2]][hk[[2]]]))

1 - sqrt( ( 1 - Gt/Ga) * (1 - Gt/Gb) )

rm(trinarized.hk)

# Remove from lists

trinarized.exp[c(1,5,7)] <- lapply(trinarized.exp[c(1,5,7)], function(y) lapply(y, function(x) x[!(names(x) %in% names(hk[[1]]))]))

trinarized.exp[c(2:4,6,8:10)] <- lapply(trinarized.exp[c(2:4,6,8:10)], function(y) lapply(y, function(x) x[!(names(x) %in% names(hk[[1]]))]))


# Save out housekeeping filtered trinarized expression

saveRDS(trinarized.exp, file = paste("trinarized_expression_hk-filtered_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))


