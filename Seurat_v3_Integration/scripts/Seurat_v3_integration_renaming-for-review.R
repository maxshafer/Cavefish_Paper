library(Seurat)
library(stringr)

# Load Seurat objects
hypo.integrated <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Hypo_integrated_127k_1500VFs_100Dims_v3.rds")

# For each object, need to change the Subtype and Subcluster names, and rename the meta data columns

head(hypo.integrated@meta.data)
#                          orig.ident nCount_RNA nFeature_RNA   sex RNA_snn_res.0.6 seurat_clusters       Subtype  SubclusterType SubclusterType_number   species RNA_snn_res.0.7 saved.idents species_Subtype
# hypo1_AAACCTGAGAGTTGGC-1      hypo1       1190          563 male1               2               4 Progenitors_1 Progenitors_1_6                     6 zebrafish            <NA>         <NA>            <NA>
# hypo1_AAACCTGCAAGCGTAG-1      hypo1       1979          882 male1               2               4 Progenitors_1 Progenitors_1_3                     3 zebrafish            <NA>         <NA>            <NA>
# hypo1_AAACCTGCAGTAGAGC-1      hypo1       1954         1025 male1              13               2        Glut_5        Glut_5_0                     0 zebrafish            <NA>         <NA>            <NA>
# hypo1_AAACCTGCATCACCCT-1      hypo1       1662          813 male1               2               4 Progenitors_1 Progenitors_1_5                     5 zebrafish            <NA>         <NA>            <NA>
# hypo1_AAACCTGGTAAGAGAG-1      hypo1       1240          575 male1              14              16   Endothelial   Endothelial_0                     0 zebrafish            <NA>         <NA>            <NA>
# hypo1_AAACCTGGTAGAGCTG-1      hypo1       1626          921 male1               9               0        GABA_2        GABA_2_0                     0 zebrafish            <NA>         <NA>            <NA>
#                          species_SubclusterType morph_Subtype morph_SubclusterType integrated_snn_res.0.1 integrated_snn_res.0.15 integrated_snn_res.0.2 integrated_snn_res.0.3 integrated_Subtype species.2
# hypo1_AAACCTGAGAGTTGGC-1                   <NA>          <NA>                 <NA>                      2                       3                      2                      4        Progenitors zebrafish
# hypo1_AAACCTGCAAGCGTAG-1                   <NA>          <NA>                 <NA>                      2                       3                      2                      4        Progenitors zebrafish
# hypo1_AAACCTGCAGTAGAGC-1                   <NA>          <NA>                 <NA>                      0                       1                      3                      2             Glut_1 zebrafish
# hypo1_AAACCTGCATCACCCT-1                   <NA>          <NA>                 <NA>                      2                       3                      2                      4        Progenitors zebrafish
# hypo1_AAACCTGGTAAGAGAG-1                   <NA>          <NA>                 <NA>                     11                      14                     16                     16        Endothelial zebrafish
# hypo1_AAACCTGGTAGAGCTG-1                   <NA>          <NA>                 <NA>                      0                       0                      0                      0             GABA_0 zebrafish
#                          integrated_SubclusterType integrated_SubclusterType_number     morph integrated_Subtype_species integrated_SubclusterType_species percent.mt
# hypo1_AAACCTGAGAGTTGGC-1             Progenitors_1                                1 zebrafish      Progenitors_zebrafish           Progenitors_1_zebrafish  7.8151261
# hypo1_AAACCTGCAAGCGTAG-1            Progenitors_10                               10 zebrafish      Progenitors_zebrafish          Progenitors_10_zebrafish  0.5053057
# hypo1_AAACCTGCAGTAGAGC-1                  Glut_1_1                                1 zebrafish           Glut_1_zebrafish                Glut_1_1_zebrafish  2.7635619
# hypo1_AAACCTGCATCACCCT-1             Progenitors_5                                5 zebrafish      Progenitors_zebrafish           Progenitors_5_zebrafish  1.2033694
# hypo1_AAACCTGGTAAGAGAG-1             Endothelial_0                                0 zebrafish      Endothelial_zebrafish           Endothelial_0_zebrafish  3.3870968
# hypo1_AAACCTGGTAGAGCTG-1                  GABA_0_1                                1 zebrafish           GABA_0_zebrafish                GABA_0_1_zebrafish  4.2435424


Idents(hypo.integrated) <- "integrated_Subtype"

# This is the current cluster IDs arranged by cluster size (for renaming to "Neuronal_")

current.cluster.ids <- c("GABA_0","GABA_1","Progenitors","Prdx1_Positive",
                         "Glut_0","Glut_1","Glut_2","GABA_2","Leucocytes",     
                         "Glut_3","Microglia","Glut_4","Glut_5","Macrophages",
                         "GABA_3","Oligodendrocytes","OPCs","Endothelial",     
                         "GABA_4","Glut_6","Erythrocytes","Cilliated","Lymphatic",
                         "Ependymal","Glut_7" )


new.cluster.ids <- c("Neuronal_00","Neuronal_01","Progenitors","Neuronal_02",
                     "Neuronal_03","Neuronal_04","Neuronal_05","Neuronal_06","Leucocytes",     
                     "Neuronal_07","Microglia","Neuronal_08","Neuronal_09","Macrophages",
                     "Neuronal_10","Oligodendrocytes","OPCs","Endothelial",     
                     "Neuronal_11","Neuronal_12","Erythrocytes","Cilliated","Lymphatic",
                     "Ependymal","Neuronal_13" )

Idents(hypo.integrated) <- factor(new.cluster.ids[match(Idents(hypo.integrated), current.cluster.ids)], levels = c("Endothelial","Erythrocytes","Ciliated","Ependymal","Progenitors","OPCs","Oligodendrocytes", sort(new.cluster.ids[grep("Neuronal_", new.cluster.ids)]), "Lymphatic","Leucocytes","Macrophages","Microglia"))
hypo.integrated$integrated_Cluster <- Idents(hypo.integrated)

hypo.integrated@meta.data$integrated_Subcluster <- paste(hypo.integrated@meta.data$integrated_Cluster, hypo.integrated@meta.data$integrated_SubclusterType_number, sep = "_")

## Set factor levels for Subclusters based on Cluster factor levels

index <- as.data.frame(unique(hypo.integrated@meta.data[,c("integrated_Cluster", "integrated_Subcluster", "integrated_SubclusterType_number")]))
index <- index[order(index[,1], index[,3]),]

hypo.integrated@meta.data$integrated_Subcluster <- factor(hypo.integrated@meta.data$integrated_Subcluster, levels = index$integrated_Subcluster)

saveRDS(hypo.integrated, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Hypo_integrated_127k_1500VFs_100Dims_vR.rds")


