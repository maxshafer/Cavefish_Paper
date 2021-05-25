library(Seurat)
library(stringr)

# Load Seurat objects
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")

# For each object, need to change the Subtype and Subcluster names, and rename the meta data columns

head(hypo.ast@meta.data)
#                         orig.ident nCount_RNA nFeature_RNA          species  sex RNA_snn_res.0.6 seurat_clusters RNA_snn_res.0.7  saved.idents
# SFM2_AAACCTGAGAATGTTG-1       SFM2       3796         1173 astyanax_surface male               3               3               3 Progenitors_1
# SFM2_AAACCTGCAATGGAAT-1       SFM2       1073          539 astyanax_surface male               3               3               3 Progenitors_1
# SFM2_AAACCTGCAATGTTGC-1       SFM2        862          559 astyanax_surface male              22              23              23        GABA_5
# SFM2_AAACCTGCAGGTTTCA-1       SFM2        731          436 astyanax_surface male               0               0               0        Glut_0
# SFM2_AAACCTGGTACGCACC-1       SFM2        961          438 astyanax_surface male               4               4               4     Microglia
# SFM2_AAACCTGGTATCTGCA-1       SFM2       1543          753 astyanax_surface male               8               7               7        Glut_2
#                               Subtype   SubclusterType SubclusterType_number                species_Subtype            species_SubclusterType
# SFM2_AAACCTGAGAATGTTG-1 Progenitors_1 Progenitors_1_10                    10 Progenitors_1_astyanax_surface Progenitors_1_10_astyanax_surface
# SFM2_AAACCTGCAATGGAAT-1 Progenitors_1  Progenitors_1_0                     0 Progenitors_1_astyanax_surface  Progenitors_1_0_astyanax_surface
# SFM2_AAACCTGCAATGTTGC-1        GABA_5         GABA_5_0                     0        GABA_5_astyanax_surface         GABA_5_0_astyanax_surface
# SFM2_AAACCTGCAGGTTTCA-1        Glut_0         Glut_0_0                     0        Glut_0_astyanax_surface         Glut_0_0_astyanax_surface
# SFM2_AAACCTGGTACGCACC-1     Microglia      Microglia_0                     0     Microglia_astyanax_surface      Microglia_0_astyanax_surface
# SFM2_AAACCTGGTATCTGCA-1        Glut_3         Glut_3_2                     2        Glut_3_astyanax_surface         Glut_3_2_astyanax_surface
#                                      morph_Subtype          morph_SubclusterType        morph       Cluster       Subcluster
# SFM2_AAACCTGAGAATGTTG-1 Progenitors_1_Choy_surface Progenitors_1_10_Choy_surface Choy_surface Progenitors_1 Progenitors_1_10
# SFM2_AAACCTGCAATGGAAT-1 Progenitors_1_Choy_surface  Progenitors_1_0_Choy_surface Choy_surface Progenitors_1  Progenitors_1_0
# SFM2_AAACCTGCAATGTTGC-1        GABA_5_Choy_surface         GABA_5_0_Choy_surface Choy_surface   Neuronal_14    Neuronal_14_0
# SFM2_AAACCTGCAGGTTTCA-1        Glut_0_Choy_surface         Glut_0_0_Choy_surface Choy_surface   Neuronal_00    Neuronal_00_0
# SFM2_AAACCTGGTACGCACC-1     Microglia_Choy_surface      Microglia_0_Choy_surface Choy_surface     Microglia      Microglia_0
# SFM2_AAACCTGGTATCTGCA-1        Glut_3_Choy_surface         Glut_3_2_Choy_surface Choy_surface   Neuronal_05    Neuronal_05_2


Idents(hypo.ast) <- "Subtype"

# This is the current cluster IDs arranged by cluster size (for renaming to "Neuronal_")

current.cluster.ids <- c("Glut_0","GABA_0","GABA_Prdx1_1","Progenitors_1",
                         "Microglia","Tcells", "Glut_1", "Glut_2", "Glut_3", 
                         "Macrophages", "Bcells", "Glut_4", "Oligodendrocytes_1", 
                         "Galanin", "Otpa/b_1", "GABA_1", "GABA_2", "GABA_3", 
                         "Erythrocytes", "Ciliated/Ventricular", "Otpa/b_2", 
                         "GABA_4", "GABA_5", "OPCs", "GABA_6", "GABA_7", "Progenitors_2", 
                         "Mast_cells", "Ependymal", "Lymphatic", "GABA_8", "Endothelial", 
                         "GABA_9", "Otpa/b_3", "Neutrophils", "Oligodendrocytes_2")


new.cluster.ids <- c("Neuronal_00","Neuronal_01","Neuronal_02","Progenitors_1",
                     "Microglia","Tcells", "Neuronal_03", "Neuronal_04", "Neuronal_05", 
                     "Macrophages", "Bcells", "Neuronal_06", "Oligodendrocytes_1", 
                     "Neuronal_07", "Neuronal_08", "Neuronal_09", "Neuronal_10", "Neuronal_11", 
                     "Erythrocytes", "Ciliated", "Neuronal_12", 
                     "Neuronal_13", "Neuronal_14", "OPCs", "Neuronal_15", "Neuronal_16", "Progenitors_2", 
                     "Mast_cells", "Ependyma_cells", "Lymphatic", "Neuronal_17", "Endothelial", 
                     "Neuronal_18", "Neuronal_19", "Neutrophils", "Oligodendrocytes_2")

Idents(hypo.ast) <- factor(new.cluster.ids[match(Idents(hypo.ast), current.cluster.ids)], levels = c("Endothelial","Erythrocytes","Ciliated","Ependyma_cells","Progenitors_1","Progenitors_2","OPCs","Oligodendrocytes_1", "Oligodendrocytes_2", sort(new.cluster.ids[grep("Neuronal_", new.cluster.ids)]), "Lymphatic","Tcells", "Bcells", "Mast_cells","Neutrophils","Macrophages","Microglia"))
hypo.ast$Cluster <- Idents(hypo.ast)

hypo.ast@meta.data$Subcluster <- paste(hypo.ast@meta.data$Cluster, hypo.ast@meta.data$SubclusterType_number, sep = "_")

## Set factor levels for Subclusters based on Cluster factor levels

index <- as.data.frame(unique(hypo.ast@meta.data[,c("Cluster", "Subcluster", "SubclusterType_number")]))
index <- index[order(index[,1], index[,3]),]

hypo.ast@meta.data$Subcluster <- factor(hypo.ast@meta.data$Subcluster, levels = index$Subcluster)

# Remove and rename other columns
hypo.ast@meta.data <- hypo.ast@meta.data[,c(1:8,12:19)]
colnames(hypo.ast@meta.data) <- c("orig.ident","nCount_RNA","nFeature_RNA","species","sex","RNA_snn_res.0.6","seurat_clusters",
                                  "RNA_snn_res.0.7","Subcluster_number","species_Cluster","species_Subcluster","morph_Cluster",         
                                  "morph_Subcluster","morph","Cluster","Subcluster")

# Rename cluster and subcluster ids by morph/species
hypo.ast@meta.data$morph_Cluster <- paste(hypo.ast@meta.data$Cluster, hypo.ast@meta.data$morph, sep = "_")
hypo.ast@meta.data$species_Cluster <- paste(hypo.ast@meta.data$Cluster, hypo.ast@meta.data$species, sep = "_")
hypo.ast@meta.data$morph_Subcluster <- paste(hypo.ast@meta.data$Subcluster , hypo.ast@meta.data$morph, sep = "_")
hypo.ast@meta.data$species_Subcluster <- paste(hypo.ast@meta.data$Subcluster , hypo.ast@meta.data$species, sep = "_")

# Calcualte UMAP embedding
hypo.ast <- RunUMAP(object = hypo.ast, reduction = "pca", dims = 1:50, reduction.name = "umap", reduction.key = "umap_", seed.use = 1, check_duplicates = F, min.dist = 0.5)

saveRDS(hypo.ast, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")


