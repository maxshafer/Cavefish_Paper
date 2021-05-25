library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)
library(viridis)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

# Load and rename the subsets? also the meta data??
# Need to also set factor levels for the idents, so they dot plots look good
subsets.all <- readRDS("Hypo_integrated_128k_1500VFs_100Dims_subsets_dims1.rds")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")

# Rename subsets, and change order
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
                     "Neuronal_11","Neuronal_12","Erythrocytes","Ciliated","Lymphatic",
                     "Ependymal","Neuronal_13" )

names(subsets.all) <- new.cluster.ids[match(names(subsets.all), current.cluster.ids)]
subsets.all <- subsets.all[levels(hypo.integrated@meta.data$integrated_Cluster)]

## Rename meta data
cluster <- levels(hypo.integrated@meta.data$integrated_Cluster)
subcluster <- levels(hypo.integrated@meta.data$integrated_Subcluster)

for (i in 1:length(subsets.all)) {
  # remove cells from weird integration clusters
  subsets.all[[i]] <- subset(subsets.all[[i]], cells = WhichCells(subsets.all[[i]])[WhichCells(subsets.all[[i]]) %in% WhichCells(hypo.integrated)])
  
  subsets.all[[i]]@meta.data$integrated_Cluster <- hypo.integrated@meta.data$integrated_Cluster[match(rownames(subsets.all[[i]]@meta.data), rownames(hypo.integrated@meta.data))]
  subsets.all[[i]]@meta.data$integrated_Cluster <- factor(subsets.all[[i]]@meta.data$integrated_Cluster, levels = cluster[cluster %in% unique(subsets.all[[i]]@meta.data$integrated_Cluster)])
  
  subsets.all[[i]]@meta.data$integrated_Subcluster <- hypo.integrated@meta.data$integrated_Subcluster[match(rownames(subsets.all[[i]]@meta.data), rownames(hypo.integrated@meta.data))]
  subsets.all[[i]]@meta.data$integrated_Subcluster <- factor(subsets.all[[i]]@meta.data$integrated_Subcluster, levels = subcluster[subcluster %in% unique(subsets.all[[i]]@meta.data$integrated_Subcluster)])
}


for (i in 1:length(subsets.all)) {
  Idents(subsets.all[[i]]) <- "integrated_Subcluster"
}

markers <- list()
for(i in 1:length(subsets.all)) {
  markers[[i]] <- FindAllMarkers(subsets.all[[i]])
}

gene.lists.pos <- readRDS("drift_gene_lists_pos_trinarized_a1.5_b2_f0.1.rds")

markers.2 <- list()
plots <- list()
for(i in 1:length(subsets.all)) {
  markers.2[[i]] <- markers[[i]] %>% group_by(cluster) %>% top_n(2, avg_logFC)
  markers.2[[i]] <- markers.2[[i]]$gene
  plots[[i]] <- DotPlot(subsets.all[[i]], features = unique(markers.2[[i]]), group.by = "integrated_Subcluster", scale.max = 200, scale.min = 25) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 45, hjust = 1, vjust = 0)) + scale_colour_viridis(limits = c(-2.5,2.5))
  plots[[i]] <- plots[[i]] + theme(axis.title = element_blank(), axis.text = element_text(size = 6))
}

# Exchange hcrt (which is 3rd best) with pmchl (which is 2nd) for Neuronal_01_10
markers.2[[9]][22] <- "hcrt"
plots[[9]] <- DotPlot(subsets.all[[9]], features = unique(markers.2[[9]]), group.by = "integrated_Subcluster", scale.max = 200, scale.min = 25) + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text.y = element_text(angle = 45, hjust = 1, vjust = 0)) + scale_colour_viridis(limits = c(-2.5,2.5))
plots[[9]] <- plots[[9]] + theme(axis.title = element_blank(), axis.text = element_text(size = 6))

# This could work, say 12 and 13, or something similar?
# 2 cols x 4 rows, angle the x.axis text
# Big ones: 25 (microglia), 23 (leucs), 5 (progs), Neuronal 0, 1, 2, 3, ~ 4, 6, (could also do 5 + 7)

pdf("Figures/Hypo_integrated_dotplots_Subcluster_markers_1.pdf", height = 10.0787, width = 8.69291) 
wrap_plots(plots[c(25, 23, 5, 8:14)]) + plot_layout(guides = "collect", nrow = 5)
dev.off()

pdf("Figures/Hypo_integrated_dotplots_Subcluster_markers_2.pdf", height = 10.0787, width = 8.69291) 
wrap_plots(plots[c(15:22,24,1:4,6:7)]) + plot_layout(guides = "collect", nrow = 5)
dev.off()


