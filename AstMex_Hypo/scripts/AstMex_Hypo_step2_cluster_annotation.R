library(Seurat)
library(Matrix)
library(plyr)
library(stringr)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

## Load objects
load("AstMex_63k.Robj")
hypo <- hypo.ast

# Genes used for calling clusters (from marker gene analysis of previous clustering)
genes <- c("wu:fj39g12", "pvalb7", "fkbp1ab", "chgb", "her15.1", "zgc:165461", "prdx1", "calb2a", "pitx2", "pfn1", "cd74a", "npc2", "apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b", "cd28l", "sla2", 
           "rtn4rl2a", "meis2a", "IGKC", "plp1b", "mpz", "zic2b", "zic4", "cbln1", "adcyap1b", "galn", "LHX1 (1 of many).1", "chodl", "rergla", "otpb", "ebf3a", "uncx", "penkb", "vip", "vim", 
           "arl13b", "pde6d", "tac3a", "trh", "hbaa2", "ENSAMXG00000016444", "sst1.1", "cort", "NPY", "lhx9", "jakmip1", "nr4a2a", "ppp1r14ab", "phactr3b", "timp4.1", "tspan18a", "sfrp5", 
           "ENSAMXG00000005497", "gfap", "srgn", "rel", "p2ry11", "epd", "fetub", "mdka", "mrc1a", "stab1", "mb", "abcb4", "hopx", "slc6a3", "th", "pomca", "cldn7b", "oxt", "thnsl2", "cotl1", 
           "pgd", "ENSAMXG00000013843", "ENSAMXG00000006950", "ltb4r", "alox5ap")

DotPlot(hypo, features = genes) + RotatedAxis()
DotPlot(hypo, features = c("gng3", "slc17a6a", "gad1b", "galn", "oxt", "cd74b")) + RotatedAxis()

Idents(hypo) <- "RNA_snn_res.0.6"
current.cluster.ids <- c(0:(length(levels(Idents(hypo)))-1))

# Combining cluster 6 and 24 into GABA_1 (from marker genes)
new.cluster.ids <- c("Glut_0", "GABA_0", "GABA_Prdx1_1", "Progenitors_1",  
                     "Microglia", "Tcells", "GABA_1", "Glut_1", "Glut_2",
                     "Macrophages", "Bcells", "Glut_3", "Oligodendrocytes_1",
                     "Galanin", "Otpa/b_1", "GABA_2", "GABA_3", "Glut_4", 
                     "Erythrocytes", "Ciliated/Ventricular", "Otpa/b_2", "GABA_4", 
                     "GABA_5", "OPCs", "GABA_1", "GABA_6", "Progenitors_2", "Mast_cells",
                     "Lymphatic", "GABA_7", "Endothelial", "Otpa/b_3", 
                     "CONT", "Thrombocytes", "Neutrophils", "Oligodendrocytes_2")

Idents(hypo) <- new.cluster.ids[match(Idents(hypo), current.cluster.ids)]

hypo$Subtype <- Idents(hypo)

hypo.ast.63.2 <- hypo

saveRDS(hypo.ast.63.2, file = "AstMex_63.2k.rds")

## Remove contaminating clusters and save as different object
## One whole subtype, plus one subcluster

Idents(hypo) <- "Subtype"

# Remove all 3 CONT, or keep everything else

hypo <- subset(hypo, idents = levels(Idents(hypo))[!grepl("CONT", levels(Idents(hypo)))])

hypo <- SetAllIdent(hypo, id = "SubclusterType")

# Set factor levels for Subtypes

hypo@meta.data$Subtype <- factor(hypo@meta.data$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated/Ventricular", "Ependymal", "Progenitors_1", "Progenitors_2", 
                                                                    "OPCs", "Oligodendrocytes_1", "Oligodendrocytes_2", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", 
                                                                    "GABA_5", "GABA_6", "GABA_7", "GABA_Prdx1_1", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Galanin", 
                                                                    "Otpa/b_1", "Otpa/b_2", "Otpa/b_3", "Lymphatic", "Tcells", "Bcells", "Mast_cells", "Thrombocytes", "Neutrophils", 
                                                                    "Macrophages", "Microglia"))

# Save as new, reduced file
hypo.ast.63 <- hypo
saveRDS(hypo.ast.63, file = "AstMex_63k.rds")
