library(Matrix)
library(dplyr)
library(ggplot2)
library(cowplot)
library(tidyr)
library(RDAVIDWebService)
library(org.Dr.eg.db)
library(reshape2)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(pubr)
library(patchwork)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj")

load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_64k.Robj")



# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?

prop.table <- as.data.frame(table(hypo.zeb@meta.data$Type2, hypo.zeb@meta.data$species))
prop.table$ClusterNumb <- apply(prop.table, 1, function(x) length(unique(hypo.zeb@meta.data[hypo.zeb@meta.data$Type2 == x[1], "integrated_Subtype"])))
prop.table$SubclusterNumb <- apply(prop.table, 1, function(x) length(unique(hypo.zeb@meta.data[hypo.zeb@meta.data$Type2 == x[1], "integrated_SubclusterType"])))
prop.table$Proportion <- (prop.table$Freq/sum(prop.table$Freq))*100
prop.table.zeb <- prop.table

prop.table <- as.data.frame(table(hypo.ast@meta.data$Type2, hypo.ast@meta.data$species.2))
prop.table$ClusterNumb <- apply(prop.table, 1, function(x) length(unique(hypo.ast@meta.data[hypo.ast@meta.data$Type2 == x[1], "integrated_Subtype"])))
prop.table$SubclusterNumb <- apply(prop.table, 1, function(x) length(unique(hypo.ast@meta.data[hypo.ast@meta.data$Type2 == x[1], "integrated_SubclusterType"])))
prop.table$Proportion <- (prop.table$Freq/sum(prop.table$Freq))*100
prop.table.ast <- prop.table

subcluster.table <- rbind(prop.table.zeb, prop.table.ast)

subcluster.table$Var1 <- factor(subcluster.table$Var1, levels = c("Microglia", "Macrophages", "Leucocytes", "Blood", "Lymphatic", "Endothelial", "Ependymal", "Oligodendrocytes", "Oligodendrocyte_Precursor_Cells", "Progenitors", "Prdx1+_GABAergic", "Otpa/b+_neurons", "Galanin+_neurons", "Glutaminergic", "GABAergic"))

 


subcluster.plot <- ggplot(subcluster.table, aes(x = Var1, size = Proportion, y = SubclusterNumb, color = Var2)) + geom_point() + coord_flip()

subcluster.plot + theme(axis.title = element_blank(), axis.text = element_text(size = 8), text = element_text(size = 8))

## 

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

# Make tables of cell type proportions, use 90% cutoff for Ast/Zeb specific?

prop.table <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species.2)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$high <- apply(prop.table, 1, function(x) max(x[1], x[2]))
prop.table$SubclusterType <- row.names(prop.table)
prop.table$Subtype <- hypo.integrated@meta.data$integrated_Subtype[match(prop.table$SubclusterType, hypo.integrated@meta.data$integrated_SubclusterType)]
prop.table$Neuronal <- hypo.integrated@meta.data$neuronal[match(prop.table$SubclusterType, hypo.integrated@meta.data$integrated_SubclusterType)]
prop.table$Type2 <- hypo.ast@meta.data$Type2[match(prop.table$SubclusterType, hypo.ast@meta.data$integrated_SubclusterType)]

prop.table$Type2 <- factor(prop.table$Type2, levels = c("Microglia", "Macrophages", "Leucocytes", "Blood", "Lymphatic", "Endothelial", "Ependymal", "Oligodendrocytes", "Oligodendrocyte_Precursor_Cells", "Progenitors", "Prdx1+_GABAergic", "Otpa/b+_neurons", "Galanin+_neurons", "Glutaminergic", "GABAergic"))
ggplot(prop.table, aes(x = Type2, y = high)) + geom_boxplot() + geom_jitter() + coord_flip()
