library(Seurat)
library(stringr)

# Load Seurat objects
hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")

# For each object, need to change the Subtype and Subcluster names, and rename the meta data columns

head(hypo.zeb@meta.data)
# orig.ident nCount_RNA nFeature_RNA   sex RNA_snn_res.0.6 seurat_clusters       Subtype  SubclusterType SubclusterType_number   species percent.mt
# hypo1_AAACCTGAGAGTTGGC-1      hypo1       1190          563 male1               2               2 Progenitors_1 Progenitors_1_6                     6 zebrafish  7.8151261
# hypo1_AAACCTGCAAGCGTAG-1      hypo1       1979          882 male1               2               2 Progenitors_1 Progenitors_1_3                     3 zebrafish  0.5053057
# hypo1_AAACCTGCAGTAGAGC-1      hypo1       1954         1025 male1              13              13        Glut_5        Glut_5_0                     0 zebrafish  2.7635619
# hypo1_AAACCTGCATCACCCT-1      hypo1       1662          813 male1               2               2 Progenitors_1 Progenitors_1_5                     5 zebrafish  1.2033694
# hypo1_AAACCTGGTAAGAGAG-1      hypo1       1240          575 male1              14              14   Endothelial   Endothelial_0                     0 zebrafish  3.3870968
# hypo1_AAACCTGGTAGAGCTG-1      hypo1       1626          921 male1               9               9        GABA_2        GABA_2_0                     0 zebrafish  4.2435424

Idents(hypo.zeb) <- "Subtype"

# This is the current cluster IDs arranged by cluster size (for renaming to "Neuronal_")

current.cluster.ids <- c("GABA_0", "Glut_0", "Progenitors_1", "GABA_Prdx1_1",  
                     "Glut_1", "Otpa/b_1", "Galanin", "GABA_1", "Glut_2",
                     "GABA_2", "Glut_3", "Glut_4",
                     "GABA_3", "Glut_5", "Endothelial", "GABA_Prdx1_2", "OPCs", 
                     "Microglia", "Erythrocytes", "GABA_4", "GABA_5", 
                     "Tcells", "Oligodendrocytes", "Lymphatic", "Macrophages", "GABA_6",
                     "Ependymal", "Progenitors_2", "Progenitors_3", 
                     "Glut_6", "Glut_7", "GABA_Prdx1_3", "GABA_7",
                     "vSMC/Pericytes", "GABA_8", "GABA_9")

new.cluster.ids <- c("Neuronal_00", "Neuronal_01", "Progenitors_1", "Neuronal_02",  
                     "Neuronal_03", "Neuronal_04", "Neuronal_05", "Neuronal_06", "Neuronal_07",
                     "Neuronal_08", "Neuronal_09", "Neuronal_10",
                     "Neuronal_11", "Neuronal_12", "Endothelial", "Neuronal_13", "OPCs", 
                     "Microglia", "Erythrocytes", "Neuronal_14", "Neuronal_15", 
                     "Tcells", "Oligodendrocytes", "Lymphatic", "Macrophages", "Neuronal_16",
                     "Ependymal", "Progenitors_2", "Progenitors_3", 
                     "Neuronal_17", "Neuronal_18", "Neuronal_19", "Neuronal_20",
                     "vSMC/Pericytes", "Neuronal_21", "Neuronal_22")

Idents(hypo.zeb) <- factor(new.cluster.ids[match(Idents(hypo.zeb), current.cluster.ids)], levels = c("Endothelial","vSMC/Pericytes","Erythrocytes","Ependymal","Progenitors_1","Progenitors_2","Progenitors_3","OPCs","Oligodendrocytes", sort(new.cluster.ids[grep("Neuronal_", new.cluster.ids)]), "Lymphatic","Tcells","Macrophages","Microglia"))
hypo.zeb$Cluster <- Idents(hypo.zeb)

hypo.zeb@meta.data$Subcluster <- paste(hypo.zeb@meta.data$Cluster, hypo.zeb@meta.data$SubclusterType_number, sep = "_")

## Set factor levels for Subclusters based on Cluster factor levels

index <- as.data.frame(unique(hypo.zeb@meta.data[,c("Cluster", "Subcluster", "SubclusterType_number")]))
index <- index[order(index[,1], index[,3]),]

hypo.zeb@meta.data$Subcluster <- factor(hypo.zeb@meta.data$Subcluster, levels = index$Subcluster)

# Remove and rename other columns
hypo.zeb@meta.data <- hypo.zeb@meta.data[,c(1:6,9:13)]
colnames(hypo.zeb@meta.data) <- c("orig.ident","nCount_RNA","nFeature_RNA","sex","RNA_snn_res.0.6","seurat_clusters","Subcluster_number",
                                  "species","percent.mt","Cluster","Subcluster")

# Calcualte UMAP embedding
hypo.zeb <- RunUMAP(object = hypo.zeb, reduction = "pca", dims = 1:50, reduction.name = "umap", reduction.key = "umap_", seed.use = 1, check_duplicates = F, min.dist = 0.5)

saveRDS(hypo.zeb, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vR.rds")


