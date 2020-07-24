library(Seurat)
library(Matrix)
library(plyr)
library(stringr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

## Load objects

hypo <- readRDS("DanRer_66k.rds")


# Genes used for calling clusters (from marker gene analysis of previous clustering)
genes <- c("gad2","slc17a6","pcp4a","calb2a","her4.2","prdx1","th2","pdyn","rasl11b","galn","hmx2","gsx1","nkx2.4a","tcf7l2","raver2","zic1","zic3","otpa","otpb","wu:fj39g12","tac1","pitx2","bhlhe22","lhx6","lhx8a","rtn4rl2a",
           "rtn4rl2b","hopx","flt1","sp8a","pax6a","traf4a","aplnrb","txn","atp1b2b","rprml","lhx9","tbr1b","bdnf","prdm12b","cart2","pfn1","CD53","cd74a","cd74b","npc2","mpeg1.1","apoeb","apoc1","hbaa1",
           "ba1","gsc","nts","htr1ab","crabp1a","plp1b","cd59","mrc1a","gng10","grn2","nebl","neurod1","epd","fxyd1","isl1l","scgn","ppp1r1c","rem1","gata3","shha","sst1.1","npy","vim",
           "gfap","olig2","kiss1","ppp1r14ab","crhb","npb","sst6","meis2b","gng8","tac3a","TPM1","aoc2","pnocb","penkb")

DotPlot(hypo, features = genes) + RotatedAxis()
DotPlot(hypo, features = c("gng3", "slc17a6b", "gad1b", "galn", "sst1.1", "cd74b", "epd")) + RotatedAxis()

Idents(hypo) <- "RNA_snn_res.0.6"
current.cluster.ids <- c(0:(length(levels(Idents(hypo)))-1))

# Combining cluster 6 and 24 into GABA_1 (from marker genes)
new.cluster.ids <- c("GABA_0", "Glut_0", "Progenitors_1", "GABA_Prdx1_1",  
                     "Glut_1", "Otpa/b_1", "Galanin", "GABA_1", "Glut_2",
                     "GABA_2", "Glut_3", "Glut_4",
                     "GABA_3", "Glut_5", "Endothelial", "GABA_Prdx1_2", "OPCs", 
                     "Microglia", "Erythrocytes", "GABA_4", "GABA_5", 
                     "Tcells", "Oligodendrocytes", "Lymphatic", "Macrophages", "CONT_Cerrebellum", "GABA_6",
                     "Ependymal", "Progenitors_2", "Progenitors_3", "CONT_Glut_1", 
                     "Glut_6", "Glut_7", "GABA_Prdx1_3", "CONT_Glut_2", "GABA_7",
                     "vSMC/Pericytes", "GABA_8", "GABA_9")

Idents(hypo) <- new.cluster.ids[match(Idents(hypo), current.cluster.ids)]

hypo$Subtype <- Idents(hypo)

DimPlot(hypo, reduction = "tsne", group.by = "Subtype", label = T) + NoLegend()

hypo.zeb.65.9 <- hypo

saveRDS(hypo.zeb.65.9, file = "DanRer_65.9k.rds")


## Remove contaminating clusters and save as different object
## One whole subtype, plus one subcluster

Idents(hypo) <- "Subtype"

# Remove all 3 CONT, or keep everything else

hypo <- subset(hypo, idents = levels(Idents(hypo))[!grepl("CONT", levels(Idents(hypo)))])

# Set factor levels for Subtypes

hypo@meta.data$Subtype <- factor(hypo@meta.data$Subtype, levels = c("Endothelial", "vSMC/Pericytes", "Erythrocytes", "Ependymal", "Progenitors_1", "Progenitors_2", "Progenitors_3", 
                                                                    "OPCs", "Oligodendrocytes", "GABA_0","GABA_1","GABA_2","GABA_3","GABA_4", "GABA_5","GABA_6","GABA_7","GABA_8", "GABA_9",
                                                                    "GABA_Prdx1_1","GABA_Prdx1_2", "GABA_Prdx1_3", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Glut_7",
                                                                    "Galanin", "Otpa/b_1", "Lymphatic", "Tcells","Macrophages", "Microglia"))

# Save as new, reduced file
hypo.zeb.65 <- hypo
saveRDS(hypo.zeb.65, file = "DanRer_65k.rds")
