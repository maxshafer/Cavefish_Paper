library(Seurat)
library(Matrix)
library(dplyr)
library(datapasta)
library(corrgram)
library(corrplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo")

load("Shafer_Hypo_66k.Robj")

# Generate list of markers for low res clustering

dir.markers <- paste("CSV/full_dataset/cluster_markers/", list.files("CSV/full_dataset/cluster_markers/"), sep = "")

markers <- vector()
for (i in 1:length(dir.markers)) {
	temp <- as.data.frame(read.csv(file=dir.markers[i], head =T))
	temp <- as.vector(temp[1:50, 8])
	markers <- c(markers,temp)
}

markers <- tibble(markers = markers)
markers <- markers %>% group_by(markers) %>% tally(sort = T)

markers <- markers[markers$n > 1,]


# Find genes which are markers for multiple subclusters

dir.markers.sub <- paste("Figures/subclustering/CSV/", list.files("Figures/subclustering/CSV/"), sep = "")

markers.sub <- vector()
for (i in 1:length(dir.markers.sub)) {
	temp <- as.data.frame(read.csv(file=dir.markers.sub[i], head =T))
	temp <- as.data.frame(temp %>% group_by(cluster) %>% top_n(-50, p_val_adj))
	temp <- as.vector(temp$gene)
	markers.sub <- c(markers.sub,temp)
}

markers.sub <- tibble(markers.sub = markers.sub)
markers.sub <- markers.sub %>% group_by(markers.sub) %>% tally(sort = T)

markers.sub <- markers.sub[markers.sub$n > 1,]

clipr::write_clip(markers.sub$markers.sub)

## Go to David, run analysis, get TFs (nuclear), use datapasta

df_paste()
tf_genes <- c("bach1b", "barhl2", "fosb", "fezf1", "fezf2","fosl1a", "h1fx", "h3f3b.1", "h3f3d", "hmx2", "isl1", "jdp2b", "lhx1a", "lhx5", "lhx6", "lhx8a", "lhx9", "lmx1bb", "meis1b", "meis2a", "nkx2.2a", "nkx2.4a", "nkx2.4b", "pou2f2a", "pou3f1", "pou3f3a", "pou3f3b", "prdm12b", "six3a", "six3b", "six6a", "six6b", "sox1a", "sox1b", "sox2", "sox3", "tbx3a", "tead1b", "tgif1", "uncx", "uncx4.1", "ascl1a", "atf3", "bsx", "dlx1a", "dlx2a", "dlx2b", "dlx5a", "dlx6a", "ebf3a", "egr1", "egr2b", "egr3", "emx1", "emx2", "emx3", "esrrga", "etv1", "etv5a", "foxa1", "foxb1a", "foxd2", "foxg1a", "foxo1a", "foxp1b", "foxp4", "fosl2", "gbx1", "histh1l", "id2b", "insm1a", "insm1b", "irx3a", "irx5a", "junba", "junbb", "jund", "jun", "mef2cb", "npas4a", "npas4b", "neurod1", "neurod6a", "neurod6b", "nfia", "nfixb", "nfil3", "nr2f1a", "nr2f2", "nr4a1", "nr4a3", "nusap1", "onecut1", "otpa", "otpb", "pax6a", "pitx2", "pbx1a", "pbx3b", "pcna", "prox1a", "rfx4", "shox2", "shox", "si:ch211-103n10.5", "si:ch211-59d15.4", "si:dkey-23a13.7", "top2a", "fosaa", "fosab", "myca", "mycb", "vax1", "zbtb16a", "zbtb18", "zfhx3", "zfhx4")

ggplot(markers.sub[1:30,], aes(x = markers.sub, y = n)) + geom_bar(stat = "identity") + coord_flip()
ggplot(markers.sub[tf_genes %in% markers.sub$markers.sub,], aes(x = markers.sub, y = n)) + geom_bar(stat = "identity")

FeaturePlot(object = SubsetData(hypo.list[[2]], ident.use = c("Otpa/b_1", "Otpa/b_2", "Otpa/b_3")), no.axes = T, reduction.use = "tsne", features.plot = c("otpb", "avp", "oxt", "trh", "th", "sst1.1", "thnsl2"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = 1)

FeaturePlot(object = subcluster, no.axes = T, reduction.use = "tsne", features.plot = c("avp"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = 3)

FeaturePlot(object = subcluster, no.axes = T, reduction.use = "tsne", features.plot = c("ENSAMXG00000025407", "ctsla"), overlay = TRUE, cols.use = c("grey85", "blue", "red", "green"), pt.size = 3, no.legend = F)

TSNEPlot(subcluster, group.by = "species_subcluster", colors.use = c("springgreen3", "springgreen4", "lightgoldenrod2", "lightgoldenrod4"), pt.size = 3, do.label = FALSE, do.return = FALSE, no.legend = FALSE)



DotPlot(object = SubsetData(hypo.list[[1]], ident.use = c("GABA_Prdx1_1")), genes.plot = c("th2", "tph1a", "prdx1"), group.by = "species_subtype", plot.legend = TRUE, x.lab.rot = TRUE, scale.min = 0, scale.max = 100, dot.scale = 3)

DotPlot(object = SubsetData(hypo.list[[2]], ident.use = c("Glut_5")), genes.plot = c("hcrt"), group.by = "species_subcluster", plot.legend = TRUE, x.lab.rot = TRUE)



# Can we look at correlation amoung these genes?

markers.sub.cut <- markers.sub[markers.sub$n > 3,]

matrix <- t(as.matrix(hypo@data[markers.sub.cut$markers.sub,]))

matrix.cor <- cor(matrix)

corrplot(matrix.cor, type = "full", method = "color", order = "hclust", tl.col = "#666666", tl.srt = 45, tl.cex = 0.5, col = brewer.pal(n=15, name="PuOr"), diag = F, addgrid.col = NA, cl.ratio = 0.07)




markers <- readRDS("Shafer_Hypo_markers.Subtype.rds")





# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
names(go_lists) <- list.files()[grep("GO", list.files())]
names(go_lists)

go_lists <- lapply(go_lists, function(x) x[x %in% row.names(hypo.zeb@data)])



DotPlot(object = hypo.zeb, genes.plot = go_lists[[1]], group.by = "Subtype", plot.legend = TRUE, x.lab.rot = TRUE)

