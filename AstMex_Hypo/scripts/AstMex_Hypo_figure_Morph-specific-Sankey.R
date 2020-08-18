library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(networkD3)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

hypo.ast <- readRDS("AstMex_63k.rds")

## New try

prop.table <- table(hypo.ast@meta.data$SubclusterType, hypo.ast@meta.data$morph)

test <- t(apply(prop.table, 1, function(x) ifelse(x/sum(x) < 0.1, 0, 1)))

shared <- test[test[,1] == 1, 2:4]

allcaves <- dim(shared[shared[,4] == 3,])
pandt <- nrow(shared[shared[,1] == 0 & shared[,2] == 1 & shared[,3] == 1,])
pandm <- nrow(shared[shared[,1] == 1 & shared[,2] == 1 & shared[,3] == 0,])
tandm <- nrow(shared[shared[,1] == 1 & shared[,2] == 0 & shared[,3] == 1,])
m <- nrow(shared[shared[,1] == 1 & shared[,2] == 0 & shared[,3] == 0,])
p <- nrow(shared[shared[,1] == 0 & shared[,2] == 1 & shared[,3] == 0,])
t <- nrow(shared[shared[,1] == 0 & shared[,2] == 0 & shared[,3] == 1,])


callcaves <- nrow(test[test[,1] == 0 & test[,2] == 1 & test[,3] == 1 & test[,4] == 1,])
cpandt <- nrow(test[test[,1] == 0 & test[,2] == 0 & test[,3] == 1 & test[,4] == 1,])
cpandm <- nrow(test[test[,1] == 0 & test[,2] == 1 & test[,3] == 1 & test[,4] == 0,])
ctandm <- nrow(test[test[,1] == 0 & test[,2] == 1 & test[,3] == 0 & test[,4] == 1,])
cm <- nrow(test[test[,1] == 0 & test[,2] == 1 & test[,3] == 0 & test[,4] == 0,])
cp <- nrow(test[test[,1] == 0 & test[,2] == 0 & test[,3] == 1 & test[,4] == 0,])
ct <- nrow(test[test[,1] == 0 & test[,2] == 0 & test[,3] == 0 & test[,4] == 1,])

choy <- nrow(test[test[,1] == 1 & test[,2] == 0 & test[,3] == 0 & test[,4] == 0,])

shared <- c(rep("Shared", 7), rep("Morph-specific", 8))
morph <- c(rep("Shared2", 7), rep("Cave-morph-specific", 7), "Surface-morph-specific")
caves <- c("All Caves", "Pachon + Tinaja", "Pachon + Molino", "Tinaja + Molino", "Pachon", "Tinaja", "Molino", "All Caves", "Pachon + Tinaja", "Pachon + Molino", "Tinaja + Molino", "Pachon", "Tinaja", "Molino", "Choy-surface-morph")
freq <- c(allcaves, pandt, pandm, tandm, p, t, m, callcaves, cpandt, cpandm, ctandm, cp, ct, cm, choy)


data <- data.frame(shared = shared, morph = morph, caves = caves, freq = freq)

### For networkD3

# Make df with Var1, Var2 and the values (var1 and var2 change)
data4 <- data[1:7,2:4]
data5 <- data[8:15,c(3,2,4)]
colnames(data4) <- c("Var.1", "Var.2", "value")
colnames(data5) <- c("Var.1", "Var.2", "value")
data4 <- rbind(data4, data5)

nodes <- unique(c(as.character(unique(data4$Var.2)), as.character(unique(data4$Var.1))))
nodes <- as.data.frame(nodes)

data4$Var1 <- as.numeric(row.names(nodes)[match(data4$Var.1, nodes$nodes)])-1
data4$Var2 <- as.numeric(row.names(nodes)[match(data4$Var.2, nodes$nodes)])-1

data4 <- data4[data4$value > 0,]

sankeyNetwork(Links = data4, Nodes = nodes, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes", units = "cells", fontSize = 8, nodeWidth = 40, nodePadding = 10, height = 800, width = 800, iterations = 10000)



