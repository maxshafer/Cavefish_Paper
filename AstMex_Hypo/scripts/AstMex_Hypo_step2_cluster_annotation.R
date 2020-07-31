library(Seurat)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

## Load objects
hypo <- readRDS("AstMex_63k.rds")

# Genes used for calling clusters (from marker gene analysis of previous clustering)
genes <- c("wu:fj39g12", "pvalb7", "fkbp1ab", "chgb", "her15.1", "zgc:165461", "prdx1", "calb2a", "pitx2", "pfn1", "cd74a", "npc2", "apoeb", "apoc1", "csf3a", "dusp2", "cxcr4b", "cd28l", "sla2", 
           "rtn4rl2a", "meis2a", "IGKC", "plp1b", "mpz", "zic2b", "zic4", "cbln1", "adcyap1b", "galn", "LHX1 (1 of many).1", "chodl", "rergla", "otpb", "ebf3a", "uncx", "penkb", "vip", "vim", 
           "arl13b", "pde6d", "tac3a", "trh", "hbaa2", "ENSAMXG00000016444", "sst1.1", "cort", "NPY", "lhx9", "jakmip1", "nr4a2a", "ppp1r14ab", "phactr3b", "timp4.1", "tspan18a", "sfrp5", 
           "ENSAMXG00000005497", "gfap", "srgn", "rel", "p2ry11", "epd", "fetub", "mdka", "mrc1a", "stab1", "mb", "abcb4", "hopx", "slc6a3", "th", "pomca", "cldn7b", "oxt", "thnsl2", "cotl1", 
           "pgd", "ENSAMXG00000013843", "ENSAMXG00000006950", "ltb4r", "alox5ap")

DotPlot(hypo, features = c("gng3", "slc17a6a", "gad1b",genes), group.by = "RNA_snn_res.0.6") + RotatedAxis()
DotPlot(hypo, features = c("gng3", "slc17a6a", "gad1b", "galn", "oxt", "cd74b")) + RotatedAxis()

Idents(hypo) <- "RNA_snn_res.0.6"

# No thrombocytes this clustering (they group with neutrophils)
hypo <- RenameIdents(hypo, '0' = "Glut_0", '1' = "GABA_0", '2' = "GABA_Prdx1_1", '3' = "Progenitors_1", '4' = "Microglia", '5' = "Tcells", '6' = "Glut_1", '7' = "Glut_2", '8' = "Glut_3", '9' = "Macrophages", 
                            '10' = "Bcells", '11' = "Glut_4", '12' = "Oligodendrocytes_1", '13' = "Galanin", '14' = "Otpa/b_1", '15' = "GABA_1", '16' = "GABA_2", '17' = "GABA_3", '18' = "Erythrocytes",
                            '19' = "Ciliated/Ventricular", '20' = "Otpa/b_2", '21' = "GABA_4", '22' = "GABA_5", '23' = "OPCs", '24' = "GABA_6", '25' = "GABA_7", '26' = "Progenitors_2", '27' = "Mast_cells", '28' = "Ependymal",
                            '29' = "Lymphatic", '30' = "GABA_8", '31' = "Endothelial", '32' = "GABA_9", '33' = "Otpa/b_3", '34' = "CONT", '35' = "Neutrophils", '36' = "Oligodendrocytes_2")

hypo$Subtype <- Idents(hypo)

hypo.ast.63.2 <- hypo

saveRDS(hypo.ast.63.2, file = "AstMex_63.2k.rds")

## Remove contaminating clusters and save as different object
## One whole subtype, plus one subcluster

Idents(hypo) <- "Subtype"

# Remove CONT cluster (from only 1 sample), and keep everything else as new smaller object

hypo <- subset(hypo, idents = levels(Idents(hypo))[!grepl("CONT", levels(Idents(hypo)))])

# Set factor levels for Subtypes

hypo@meta.data$Subtype <- factor(hypo@meta.data$Subtype, levels = c("Endothelial", "Erythrocytes", "Ciliated/Ventricular", "Ependymal", "Progenitors_1", "Progenitors_2", 
                                                                    "OPCs", "Oligodendrocytes_1", "Oligodendrocytes_2", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", 
                                                                    "GABA_5", "GABA_6", "GABA_7", "GABA_8", "GABA_9", "GABA_Prdx1_1", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Galanin", 
                                                                    "Otpa/b_1", "Otpa/b_2", "Otpa/b_3", "Lymphatic", "Tcells", "Bcells", "Mast_cells", "Neutrophils", "Macrophages", "Microglia"))

# Save as new, reduced file
hypo.ast.63 <- hypo
saveRDS(hypo.ast.63, file = "AstMex_63k.rds")
