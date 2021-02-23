library(Seurat)
library(reticulate)


# Set hyperparameters from arguments or in file
args <- commandArgs(trailingOnly = TRUE)

print(args)

a <- as.numeric(args[[1]])
b <- as.numeric(args[[2]])
f <- as.numeric(args[[3]])

# a = 1.5
# b = 2
# f = 0.1

# import scipy functions (from Zeisel et al scripts) to reproduce behaviour (particularily of incomplete beta)
scipy.special <- import("scipy.special")

# Define functions

betabinomial_trinarize <- function(raw.data = raw.data, n = n, a = 1.5, b = 2, f = 0.2) {
  n <- length(raw.data)
  k <- length(raw.data[raw.data > 0])
  
  incb <- scipy.special$betainc(a + k, b - k + n, f)
  if (incb == 0) {
    p <- 1.0
  } else {
    p <- 1.0 - exp(log(incb) + lbeta(a + k, b - k + n) + lgamma(a + b + n) - lgamma(a + k) - lgamma(b - k + n)) # scipy.special$betaln + math$lgamma if using python version, which is slower then built in functions
  }
  return(p)
}

# Load data
setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

DefaultAssay(hypo.integrated) <- "RNA"

Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

Idents(hypo.integrated.ast) <- "species"
hypo.integrated.surface <- subset(hypo.integrated.ast, idents = "astyanax_surface")
hypo.integrated.cave <- subset(hypo.integrated.ast, idents = "astyanax_cave")

# make list
trinarized.genes <- list()

############## Trinarize Clusters ##############

# Trinarize cluster gene expression

Idents(hypo.integrated) <- "integrated_Cluster"
Idents(hypo.integrated.zeb) <- "integrated_Cluster"
Idents(hypo.integrated.ast) <- "integrated_Cluster"
Idents(hypo.integrated.surface) <- "integrated_Cluster"
Idents(hypo.integrated.cave) <- "integrated_Cluster"

# In zebrafish (hypo.integrated.zeb), for each cluster find if a gene is above the trinarization threshold
# For each cluster a list of all genes and the probability? Yes, and then maybe make another data structure of only those that pass (if too big)

### For species ###

cluster.zebrafish <- lapply(levels(Idents(hypo.integrated.zeb))[!(levels(Idents(hypo.integrated.zeb)) %in% "Ciliated")], function(x) apply(GetAssayData(object = hypo.integrated.zeb, slot = "counts")[, WhichCells(hypo.integrated.zeb, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(cluster.zebrafish) <- levels(Idents(hypo.integrated.zeb))[!(levels(Idents(hypo.integrated.zeb)) %in% "Ciliated")]

cluster.astyanax <- lapply(levels(Idents(hypo.integrated.ast))[!(levels(Idents(hypo.integrated.ast)) %in% "Ciliated")], function(x) apply(GetAssayData(object = hypo.integrated.ast, slot = "counts")[, WhichCells(hypo.integrated.ast, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(cluster.astyanax) <- levels(Idents(hypo.integrated.ast))[!(levels(Idents(hypo.integrated.ast)) %in% "Ciliated")]
cluster.astyanax <- cluster.astyanax[c(names(cluster.zebrafish))]

trinarized.genes[[1]] <- cluster.zebrafish
trinarized.genes[[2]] <- cluster.astyanax

### For cave vs surface ###

cluster.surface <- lapply(levels(Idents(hypo.integrated.surface)), function(x) apply(GetAssayData(object = hypo.integrated.surface, slot = "counts")[, WhichCells(hypo.integrated.surface, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(cluster.surface) <- levels(Idents(hypo.integrated.surface))

cluster.cave <- lapply(levels(Idents(hypo.integrated.cave)), function(x) apply(GetAssayData(object = hypo.integrated.cave, slot = "counts")[, WhichCells(hypo.integrated.cave, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(cluster.cave) <- levels(Idents(hypo.integrated.cave))

trinarized.genes[[3]] <- cluster.surface
trinarized.genes[[4]] <- cluster.cave

############## Trinarize Subclusters ##############

# Trinarize cluster gene expression

Idents(hypo.integrated) <- "integrated_Subcluster"
Idents(hypo.integrated.zeb) <- "integrated_Subcluster"
Idents(hypo.integrated.ast) <- "integrated_Subcluster"

# Find species specific clusters
prop.table <- table(hypo.integrated@meta.data$integrated_Subcluster, hypo.integrated@meta.data$species.2)
#prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
zeb.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
ast.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo.integrated))
index <- subclusters[!(subclusters %in% c(zeb.names, ast.names))]

### For species ###

subcluster.zebrafish <- lapply(index, function(x) apply(GetAssayData(object = hypo.integrated.zeb, slot = "counts")[, WhichCells(hypo.integrated.zeb, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(subcluster.zebrafish) <- index

subcluster.astyanax <- lapply(index, function(x) apply(GetAssayData(object = hypo.integrated.ast, slot = "counts")[, WhichCells(hypo.integrated.ast, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(subcluster.astyanax) <- index
subcluster.astyanax <- subcluster.astyanax[c(names(subcluster.zebrafish))]

trinarized.genes[[5]] <- subcluster.zebrafish
trinarized.genes[[6]] <- subcluster.astyanax

# Find for zeb and ast specific cell types

specific.zeb <- lapply(zeb.names, function(x) apply(GetAssayData(object = hypo.integrated.zeb, slot = "counts")[, WhichCells(hypo.integrated.zeb, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(specific.zeb) <- zeb.names

specific.ast <- lapply(ast.names, function(x) apply(GetAssayData(object = hypo.integrated.ast, slot = "counts")[, WhichCells(hypo.integrated.ast, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(specific.ast) <- ast.names

trinarized.genes[[7]] <- specific.zeb
trinarized.genes[[8]] <- specific.ast

### For Cave vs Surface ###

Idents(hypo.integrated.ast) <- "integrated_Subcluster"
Idents(hypo.integrated.surface) <- "integrated_Subcluster"
Idents(hypo.integrated.cave) <- "integrated_Subcluster"

# Find morph specific clusters
prop.table <- table(hypo.integrated.ast@meta.data$integrated_Subcluster, hypo.integrated.ast@meta.data$species)
# prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
surface.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
cave.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

# Make index
subclusters <- levels(Idents(hypo.integrated.ast))
index2 <- subclusters[!(subclusters %in% c(zeb.names, surface.names, cave.names))]

subcluster.surface <- lapply(index2, function(x) apply(GetAssayData(object = hypo.integrated.surface, slot = "counts")[, WhichCells(hypo.integrated.surface, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(subcluster.surface) <- index2

subcluster.cave <- lapply(index2, function(x) apply(GetAssayData(object = hypo.integrated.cave, slot = "counts")[, WhichCells(hypo.integrated.cave, idents = x)], 1, function(y) betabinomial_trinarize(raw.data = y, a = a, b = b, f = f)))
names(subcluster.cave) <- index2
subcluster.cave <- subcluster.cave[c(names(subcluster.surface))]

trinarized.genes[[9]] <- subcluster.surface
trinarized.genes[[10]] <- subcluster.cave


## Name and save

names(trinarized.genes) <- c("cluster.zebrafish", "cluster.astyanax", "cluster.surface", "cluster.cave", "subcluster.zebrafish", "subcluster.astyanax", "specific.zebrafish", "specific.astyanax", "subcluster.surface", "subcluster.cave")

saveRDS(trinarized.genes, file = paste("trinarized_expression_a",a,"_b",b, "_f",f,".rds", sep = ""))

## Pick those that pass the expression cutoff (p > 0.95)

trinarized.exp <- lapply(trinarized.genes, function(x) lapply(x, function(y) y[y > 0.95]))

saveRDS(trinarized.exp, file = paste("trinarized_expression_a",a,"_b",b, "_f",f,"_cutoff.rds", sep = ""))

