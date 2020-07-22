pcs <- Embeddings(hypo.integrated, reduction = "pca")[,1:50]


pcs <- melt(pcs)

pcs$subcluster <-  hypo.integrated@meta.data$integrated_SubclusterType[match(pcs$Var1, row.names(hypo.integrated@meta.data))]

ggplot(pcs[pcs$subcluster == "Prdx1_Positive_11",], aes(x = Var2, y = value)) + geom_jitter() + theme(axis.text.x = element_text(angle = 90))


gene.pcs <- hypo.integrated@reductions$pca
gene.pcs <- gene.pcs[,1:50]

head(sort(abs(gene.pcs[,35]), decreasing = T), n = 20)



subsets <- readRDS(file = "/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/Hypo_integrated_130k_1500VFs_100Dims_subsets.rds")

for(i in 1:length(subsets)) {
	DefaultAssay(subsets[[i]]) <- "RNA"
}

pcs <- Embeddings(subsets[[14]], reduction = "pca")[,1:50]

pcs <- melt(pcs)

pcs$subcluster <-  hypo.integrated@meta.data$integrated_SubclusterType[match(pcs$Var1, row.names(hypo.integrated@meta.data))]

ggplot(pcs[pcs$subcluster == "Prdx1_Positive_11",], aes(x = Var2, y = value)) + geom_jitter() + theme(axis.text.x = element_text(angle = 90)) + ylim(c(-20,20))

gene.pcs <- subsets[[14]]@reductions$pca
gene.pcs <- gene.pcs[,1:50]


head(sort(abs(gene.pcs[,20]), decreasing = T), n = 20)



DotPlot(subsets[[14]], features = head(row.names(gene.lists[[1]][["Prdx1_Positive"]]), n = 20), group.by = "integrated_SubclusterType") + theme(axis.text.x = element_text(angle = 45, hjust = 1))



DotPlot(subsets[[14]], features = head(row.names(zebrafish.markers.sub[["Prdx1_Positive_2"]]), n = 20), group.by = "integrated_SubclusterType") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(subsets[[3]], features = head(row.names(zebrafish.markers.sub[["Glut_0_6"]]), n = 5), pt.size = 1)

DimPlot(subsets[[1]], group.by = "integrated_SubclusterType", pt.size = 1)
DimPlot(subsets[[1]], group.by = "species.2", pt.size = 1)



DotPlot(subsets[[14]], features = head(row.names(cavefish.markers.sub[["Prdx1_Positive_8"]]), n = 20), group.by = "integrated_SubclusterType") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

FeaturePlot(subsets[[3]], features = head(row.names(cavefish.markers.sub[["Glut_0_5"]]), n = 16), pt.size = 1)

DimPlot(subsets[[14]], group.by = "integrated_SubclusterType", pt.size = 1)


table(subsets[[4]]@meta.data$integrated_SubclusterType, subsets[[4]]@meta.data$species)



DotPlot(hypo.integrated, features = head(row.names(gene.lists[[4]][["Prdx1_Positive_3"]]), n = 20), group.by = "integrated_Subtype") + theme(axis.text.x = element_text(angle = 45, hjust = 1))

Idents(hypo.integrated) <- "species"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")

DotPlot(hypo.integrated.zeb, features = head(row.names(zebrafish.markers.sub[["Glut_0_6"]]), n = 5), group.by = "Subtype") + theme(axis.text.x = element_text(angle = 45, hjust = 1))


## subset erythrocytes and these weird clusters for tsne:
Idents(hypo.integrated) <- "integrated_SubclusterType"
weird <- subset(hypo.integrated, idents = c("Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_3", "Endothelial_4", "Endothelial_5", "Endothelial_6", "Erythrocytes_0", "Erythrocytes_1", "Erythrocytes_2", "Glut_3_7", "Glut_4_6", "Glut_4_7", "GABA_2_9", "Glut_1_4"))



gene.lists <- readRDS(file = "drift_gene_lists.rds")

names(gene.lists) <- c("conserved.markers", "zebrafish.markers", "astyanax.markers", "conserved.markers.sub", "zebrafish.markers.sub", "astyanax.markers.sub", "conserved.markers.ast", "surface.markers", "cave.markers", "conserved.markers.ast.sub", "surface.markers.sub", "cave.markers.sub")

# ## Need to remove GABA_5 markers, and rename GABA_6 to GABA_5
# ## First remove GABA_5

# for (i in 1:length(gene.lists)) {
	# gene.lists[[i]] <- gene.lists[[i]][!grepl("GABA_5", names(gene.lists[[i]]))]
# }

# for (i in 1:length(gene.lists)) {
	# names(gene.lists[[i]]) <- str_replace(names(gene.lists[[i]]), "GABA_6", "GABA_5")
# }


## Take only positive markers

gene.lists.pos <- list()
gene.lists.pos[c(1,4)] <- lapply(gene.lists[c(1,4)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$zebrafish_avg_logFC > 0 & x[[y]]$astyanax_avg_logFC > 0,]))
gene.lists.pos[c(2,3,5,6,8,9,11,12)] <- lapply(gene.lists[c(2,3,5,6,8,9,11,12)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$avg_logFC > 0,]))
gene.lists.pos[c(7,10)] <- lapply(gene.lists[c(7,10)], function(x) lapply(seq_along(x), function(y) x[[y]][x[[y]]$astyanax_surface_avg_logFC > 0 & x[[y]]$astyanax_cave_avg_logFC > 0,]))

names(gene.lists.pos) <- names(gene.lists[1:12])

for(i in 1:length(gene.lists.pos)){
	names(gene.lists.pos[[i]]) <- names(gene.lists[[i]])
}

