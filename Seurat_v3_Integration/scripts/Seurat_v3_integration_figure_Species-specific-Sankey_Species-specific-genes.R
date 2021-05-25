library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(networkD3)
library(viridis)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

## New try

prop.table <- table(hypo.integrated@meta.data$integrated_Subcluster, hypo.integrated@meta.data$species.2)
prop.table <- prop.table[prop.table[,1] > 0 | prop.table[,2] > 0,]
zeb.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[1] < 3 | x[2]/sum(x) > .9))]
ast.names <- row.names(prop.table)[apply(prop.table, 1, function(x) (x[2] < 3 | x[1]/sum(x) > .9))]

sort(zeb.names)
# [1] "Endothelial_6"    "Endothelial_7"    "Glut_1_5"         "Glut_1_7"         "Glut_2_0"         "Glut_2_2"         "Glut_2_4"         "Leucocytes_5"     "Microglia_2"      "Prdx1_Positive_2" "Prdx1_Positive_3" "Progenitors_8" 

sort(ast.names)
# [1] "Cilliated_0"        "Cilliated_1"        "Cilliated_2"        "Cilliated_3"        "Cilliated_4"        "GABA_3_5"           "Glut_0_13"          "Glut_1_6"           "Glut_2_1"           "Glut_3_0"           "Glut_3_6"           "Glut_6_5"           "Glut_7_2"           "Leucocytes_0"      
# [15] "Leucocytes_1"       "Leucocytes_11"      "Leucocytes_12"      "Leucocytes_13"      "Leucocytes_2"       "Leucocytes_4"       "Leucocytes_7"       "Leucocytes_9"       "Macrophages_2"      "Macrophages_4"      "Macrophages_5"      "Microglia_0"        "Microglia_1"        "Microglia_3"       
# [29] "Microglia_6"        "Microglia_8"        "Oligodendrocytes_2" "Oligodendrocytes_5" "Progenitors_13"


data1 <- data.frame(species = c(rep("Zebrafish-specific", 4), rep("Mexican-tetra-specific", 4)), type = c("Ciliated", "Glial", "Hematopoietic", "Neuronal", "Ciliated", "Glial", "Hematopoietic", "Neuronal"), freq = c(length(grep("Cilliated", zeb.names)), length(grep("Endothelial|Progenitor|Oligodendrocyte|Lymphatic|Ependymal", zeb.names)), length(grep("Microglia|Macrophage|Leucocyte|Erythrocyte", zeb.names)), length(grep("GABA|Prdx|Glut", zeb.names)), length(grep("Cilliated", ast.names)), length(grep("Endothelial|Progenitor|Oligodendrocyte|Lymphatic|Ependymal", ast.names)), length(grep("Microglia|Macrophage|Leucocyte|Erythrocyte", ast.names)), length(grep("GABA|Prdx|Glut", ast.names))))

data2 <- data.frame(shared = c(rep("Shared", 1), rep("Species-specific", 2)), species = c("NA", "Zebrafish-specific", "Mexican-tetra-specific"), freq = c(length(unique(hypo.integrated@meta.data$integrated_SubclusterType))-length(c(zeb.names,ast.names)), length(zeb.names), length(ast.names)))

colnames(data1) <- c("Var.1", "Var.2", "value")
colnames(data2) <- c("Var.1", "Var.2", "value")

data3 <- rbind(data1, data2)

### For networkD3

# Make df with Var1, Var2 and the values (var1 and var2 change)


nodes <- unique(c(as.character(unique(data3$Var.2)), as.character(unique(data3$Var.1))))
nodes <- as.data.frame(nodes)

data3$Var1 <- as.numeric(row.names(nodes)[match(data3$Var.1, nodes$nodes)])-1
data3$Var2 <- as.numeric(row.names(nodes)[match(data3$Var.2, nodes$nodes)])-1

data4 <- data3[data3$value > 0,]

sankeyNetwork(Links = data4, Nodes = nodes, Source = "Var1", Target = "Var2", Value = "value", NodeID = "nodes", units = "cells", fontSize = 8, nodeWidth = 40, nodePadding = 10, height = 800, width = 800, iterations = 10000, sinksRight = F)


## DotPlot for species-specific subcluster marker genes


para.dot <- DotPlot(hypo.integrated, features = c("si:ch211-241c24.4", "si:ch211-241c24.3", "zgc:172109", "zgc:172260", "zgc:172290", "hpcal1", "ENSAMXG00000014042",  "vip", "vipb", "ENSAMXG00000017498", "colec12", "ENSAMXG00000007230", "ppp3ca", "PPP3CA", "NPTX1", "ENSAMXG00000025407"), group.by = "integrated_Subcluster", dot.scale = 4) + RotatedAxis() + coord_flip() 
para.dot <- para.dot + theme(axis.text.x = element_text(size = 4, angle = 45), axis.text.y = element_text(size = 8), axis.title = element_blank()) + scale_colour_viridis()
para.dot <- para.dot + plot_layout(width = unit(c(130), c("mm")), height = unit(c(40), c("mm")))
dev.new()
para.dot

# jac2 = si:ch211-241c24.4
# jac3 = si:ch211-241c24.3
# jac6 = zgc:172109
# jac8 = zgc:172260
# jac9 = zgc:172290
# HPCAL1 = ENSAMXG00000014042
# COLLEC12 = ENSAMXG00000007230

