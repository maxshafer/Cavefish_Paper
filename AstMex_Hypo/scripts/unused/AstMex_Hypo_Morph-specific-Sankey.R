library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/")

load("AstMex_64k.Robj")

## New try

prop.table <- table(hypo.ast@meta.data$SubclusterType, hypo.ast@meta.data$species_morph)

test <- t(apply(prop.table, 1, function(x) ifelse(x/sum(x) < 0.1, 0, 1)))

shared <- test[test[,1] == 1, 2:4]

sum <- unlist(apply(shared, 1, function(x) sum(x)))
shared <- cbind(shared, sum)

dim(shared[shared[,4] == 3,])
dim(shared[shared[,4] == 2,])
dim(shared[shared[,4] == 1,])

two.shared <- shared[shared[,4] == 2,]
dim(two.shared[two.shared[,1] ==0,])
dim(two.shared[two.shared[,2] ==0,])
dim(two.shared[two.shared[,3] ==0,])


shared <- c(rep("Shared", 7), rep("Morph-specific", 8))
morph <- c(rep("Shared2", 7), rep("Cave-morph-specific", 7), "Surface-morph-specific")
caves <- c("All Caves", "Pachon + Tinaja", "Pachon + Molino", "Tinaja + Molino", "Pachon", "Tinaja", "Molino", "All Caves", "Pachon + Tinaja", "Pachon + Molino", "Tinaja + Molino", "Pachon", "Tinaja", "Molino", "Choy-surface-morph")
freq <- c(69, 29, 4, 24, 17, 1, 4, 0, 2, 0, 2, 3, 0, 5, 8)


data <- data.frame(shared = shared, morph = morph, caves = caves, freq = freq)

is_alluvia_form(as.data.frame(data), silent = TRUE)

data$shared <- factor(data$shared, levels = c("Shared", "Morph-specific"))
data$morph <- factor(data$morph, levels = c("Shared2", "Cave-morph-specific", "Surface-morph-specific"))
data$caves <- factor(data$caves, levels = c("All Caves", "Pachon + Tinaja", "Pachon + Molino", "Tinaja + Molino", "Pachon", "Tinaja", "Molino", "Choy-surface-morph"))


ggplot(as.data.frame(data), aes(y = freq, axis1 = shared, axis2 = morph, axis3 = caves)) + geom_alluvium(aes(fill = caves), width = 1/6) + geom_stratum(width = 1/6, fill = "white", color = "black") + scale_x_discrete(limits = c("Shared", "Cave-morphs", "Species-morph"), expand = c(.05, .05)) + scale_fill_brewer(type = "qual", palette = "Paired")

## Better figure, with shared cell types and morph specific cell types coming from opposite sides

data2 <- melt(data, id.vars = "freq")
data2$x <- rep(1:15, 3)
colnames(data2) <- c("freq", "x", "stratum", "alluvium")
data2$x <- factor(data2$x, levels = c("shared", "caves", "morph"))
data2$fill <- as.character(rep(c(1:7,1:8), 3))


data3 <- data2[c(1:7,23:30,31:45),]


ggplot(data3, aes(x = x, stratum = stratum, alluvium = alluvium, y = freq, label = stratum)) + stat_flow(aes(fill = fill)) + geom_stratum() + geom_text(stat = "stratum") + scale_fill_brewer(type = "qual", palette = "Paired")


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

sankeyNetwork(Links = data4, Nodes = nodes, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes", units = "cells", fontSize = 8, nodeWidth = 40, nodePadding = 10, height = 800, width = 800, iterations = 10000)



