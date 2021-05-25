library(Seurat)
library(patchwork)
library(ggplot2)
library(clustree)

# Load the initial Seurat objects

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo")

# Load the datasets from 10x

data.vec.zeb <- c("raw_data/Hypo1_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo1_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo1_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo1_4/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo2_4/outs/filtered_gene_bc_matrices/dr86/",
              "raw_data/Hypo3_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo3_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo3_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo3_4/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_1/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_2/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_3/outs/filtered_gene_bc_matrices/dr86/", 
              "raw_data/Hypo4_4/outs/filtered_gene_bc_matrices/dr86/")

names (data.vec.zeb) <- c("hypo1","hypo2","hypo3","hypo4","hypo5","hypo6","hypo7","hypo8","hypo9","hypo10","hypo11","hypo12","hypo13","hypo14","hypo15","hypo16")
input.data <- Read10X(data.vec.zeb)

# Create Seurat object, normalize, scale and find variable genes

hypo.zeb.og <- CreateSeuratObject(input.data, min.features = 200, project = "Danio_rerio")


setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo")

# Load the datasets from 10x

data.vec.ast <- c("raw_data/SF_M_2/AstMex102/", 
              "raw_data/SF_M_3/AstMex102/",
              "raw_data/SF_M_4/AstMex102/",
              "raw_data/SF_M_5/AstMex102/",
              "raw_data/SF_F_2/AstMex102/",
              "raw_data/SF_F_3/AstMex102/",
              "raw_data/SF_F_4/AstMex102/",
              "raw_data/SF_F_5/AstMex102/",
              "raw_data/MF_M_1/AstMex102/",
              "raw_data/MF_F_1/AstMex102/",
              "raw_data/TF_M_1/AstMex102/",
              "raw_data/TF_F_1/AstMex102/",
              "raw_data/PF_M_1/AstMex102/",
              "raw_data/PF_F_1/AstMex102/",
              "raw_data/PF_M_2/AstMex102/",
              "raw_data/PF_F_2/AstMex102/")

names(data.vec.ast) <- c("SFM2","SFM3","SFM4","SFM5","SFF2","SFF3","SFF4","SFF5","MFM1","MFF1","TFM1","TFF1","PFM1","PFF1","PFM2","PFF2")

input.data <- Read10X(data.vec.ast)

# Create Seurat object, normalize, scale and find variable genes

hypo.ast.og <- CreateSeuratObject(input.data, min.features = 200, project = "Astyanax_mexicanus")

# Save files for loading
saveRDS(hypo.zeb.og, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_og.rds")
saveRDS(hypo.ast.og, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_og.rds")

hypo.zeb.og <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_og.rds")
hypo.ast.og <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_og.rds")
## Run QC + make figures

hypo.zeb.og[["percent.mt"]] <- PercentageFeatureSet(hypo.zeb.og, pattern = "^mt-")
hypo.ast.og[["percent.mt"]] <- PercentageFeatureSet(hypo.ast.og, pattern = "^mt-")

## Make final plots

new_theme <- function() {
  theme(
    axis.text = element_text(size = 8),
    axis.title = element_text(size = 10),
    plot.title = element_blank(), 
    legend.position = "none"
  )
}

vln.plots <- VlnPlot(hypo.zeb.og, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of genes") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

vln.plots[[1]] <- VlnPlot(hypo.zeb.og, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of genes") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln.plots[[2]] <- VlnPlot(hypo.zeb.og, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of counts") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln.plots[[3]] <- VlnPlot(hypo.zeb.og, features = c("percent.mt"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("% mito.")

vln.plots[[4]] <- VlnPlot(hypo.ast.og, features = c("nFeature_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme() + ylab("# of genes") + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
vln.plots[[5]] <- VlnPlot(hypo.ast.og, features = c("nCount_RNA"), group.by = "orig.ident", ncol = 1, pt.size = 0) + new_theme()+ ylab("# of counts")

zeb.fs1 <- FeatureScatter(hypo.zeb.og, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.05) + new_theme() + theme(legend.position = c(0.85,0.35))
zeb.fs2 <- FeatureScatter(hypo.zeb.og, feature1 = "nCount_RNA", feature2 = "percent.mt", pt.size = 0.05) + new_theme() + theme(legend.position = c(0.75,0.55))

ast.fs1 <- FeatureScatter(hypo.ast.og, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size = 0.05) + new_theme() + theme(legend.position = c(0.85,0.35))

# Plot
layout <- "
AAAAAA
BBBBBB
CCCCCC
DDDDDD
EEEEEE
FFGGHH
FFGGHH
"
dev.new()
vln.plots + plot_layout(design = layout, guides = "collect", widths = unit(c(15), "mm"), height = unit(c(15), "mm"))

## Need to copy to affinity and rasterize one at a time (60k dots each)
layout2 <- "
FFGGHH
FFGGHH
"
dev.new()
zeb.fs1 + zeb.fs2 + ast.fs1 + plot_layout(design = layout2, guides = "collect", widths = unit(c(15), "mm"), height = unit(c(15), "mm"))


######### Make clustree plots for species-specific subclusters ############

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("Hypo_integrated_127k_1500VFs_100Dims_vR.rds")
Idents(hypo.integrated) <- "species.2"
hypo.integrated.zeb <- subset(hypo.integrated, idents = "zebrafish")
hypo.integrated.ast <- subset(hypo.integrated, idents = "astyanax")

# Load subsets, res 0.2, 0.25, and 0.4 were used for subclustering, I chose to use 0.25 in the end
# anything else might be from intial clustering, not subclustering

subsets.all <- readRDS(file = "Hypo_integrated_128k_1500VFs_100Dims_subsets_dims1.rds")

# Rename the subsets

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


# Make sure they all have the cluster resolutions (maybe some were re-done, but with only the resolution ultimately decided on)
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.15, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.2, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.25, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.3, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.35, print.output = 0))
subsets.all <- lapply(subsets.all, function(x) FindClusters(x, resolution = 0.4, print.output = 0))

subsets.all.meta <- list()
subsets.all.meta[c(1:24)] <- lapply(subsets.all[c(1:24)], function(x) x@meta.data[,c("integrated_snn_res.0.15", "integrated_snn_res.0.2", "integrated_snn_res.0.25", "integrated_snn_res.0.3", "integrated_snn_res.0.35", "integrated_snn_res.0.4")])
subsets.all.meta[c(25)] <- lapply(subsets.all[c(25)], function(x) x@meta.data[,c("RNA_snn_res.0.15", "RNA_snn_res.0.2", "RNA_snn_res.0.25", "RNA_snn_res.0.3", "RNA_snn_res.0.35", "RNA_snn_res.0.4")])

names(subsets.all.meta) <- names(subsets.all)

subsets.all.meta <- lapply(subsets.all.meta, function(x) {
  colnames(x) <- c("res_1.5", "res_2.0", "res_2.5", "res_3.0", "res_3.5", "res_4.0")
  return(x)
})

plots <- lapply(seq_along(subsets.all.meta), function(x) clustree(subsets.all.meta[[x]], prefix = "res_", edge_width = 0.5, node_size_range = c(1,6), alt_colour = "black") + ggtitle(paste(names(subsets.all.meta)[[x]])))
names(plots) <- names(subsets.all)

plots2 <- plots[c("Neuronal_02", "Neuronal_03", "Neuronal_04", "Neuronal_05", "Neuronal_07", "Neuronal_10", "Neuronal_12", "Neuronal_13")]

layout <- "
AAAAABBBBBCCCCC
DDDEEEFFFGGGHHH
"

plots3 <- wrap_plots(plots2) + plot_layout(heights = unit(c(45), c('mm')), widths = unit(c(20,20,20,12,12,12,12,12), c('mm')), design = layout, guides = "collect")

## Add in red boxes for chosen resolution + green boxes for species-specific subclusters in Affinity/Illustrator
pdf("Figures/Hypo_integrated_Clustree-plots_Species-specific.pdf", width = 14, height = 30)
plots3 & scale_edge_color_continuous(low = "black", high = "black")
dev.off()


####### Plots for UMI / mito / counts for species-specific cell types ########

## Maybe make prop plots for the contribution of each sample to the species-specific clusters
## Only two come from a single Mexican tetra sample, and 1 of those is the likely spatial difference one
## Can also say that it is impossible to know if they are morph/sex specific vs dissection issue (given single morph per sex)

sp.sp.subclusters <- c(c(2,3),13,c(5,6,7),c(0,1),c(0,5),5,5,2)
sp <- c(c("zebrafish", "zebrafish"), "tetra", c("zebrafish", "tetra", "zebrafish"), c("zebrafish", "tetra"), c("tetra", "tetra"), "tetra", "tetra", "tetra")

tables <- lapply(subsets.all, function(x) as.data.frame(table(x@meta.data$orig.ident, x@meta.data$integrated_snn_res.0.25)))
tables <- tables[c("Neuronal_02", "Neuronal_02", "Neuronal_03", "Neuronal_04", "Neuronal_04", "Neuronal_04", "Neuronal_05", "Neuronal_05", "Neuronal_07", "Neuronal_07", "Neuronal_10", "Neuronal_12", "Neuronal_13")]
tables <- lapply(seq_along(tables), function(x) {
  return <- tables[[x]][tables[[x]]$Var2 == sp.sp.subclusters[[x]],]
  if(sp[[x]] == "zebrafish") {
    return <- return[grep("hypo",return$Var1),]
  } else {
    return <- return[grep("F",return$Var1),]
  }
  return(return)
  })

names(tables) <- c("Neuronal_02", "Neuronal_02", "Neuronal_03", "Neuronal_04", "Neuronal_04", "Neuronal_04", "Neuronal_05", "Neuronal_05", "Neuronal_07", "Neuronal_07", "Neuronal_10", "Neuronal_12", "Neuronal_13")

prop.plots <- list()
prop.plots[c(1)] <- lapply(c(1), function(x) ggplot(tables[[x]], aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat="identity") + theme_classic() + theme(axis.title.x = element_text(size = 10)) + scale_fill_manual(values = cols3) + guides(fill = guide_legend(title = "Sample ID")) + ylab("Sample Cluster frequency") + xlab(paste(names(tables)[[x]])))
prop.plots[c(2:13)] <- lapply(c(2:13), function(x) ggplot(tables[[x]], aes(x=Var2, y=Freq, fill=Var1)) + geom_bar(stat="identity") + theme_classic() + theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 10)) + scale_fill_manual(values = cols3) + guides(fill = guide_legend(title = "Sample ID")) + ylab("Sample Cluster frequency") + xlab(paste(names(tables)[[x]])))

prop.plot <- wrap_plots(prop.plots[c(grep("zebrafish", sp), grep("tetra", sp))]) + plot_layout(nrow = 1, guides = "collect", heights = unit(c(30), c('mm')), widths = unit(c(10), c('mm')))

pdf("Figures/Hypo_integrated_Prop-plots_Species-specific.pdf", width = 14, height = 15)
prop.plot
dev.off()

## Plus plot UMIs/genes per cluster for those clusters to show they aren't abberant
## Add dataset means as dashed line?
## The zeb and tetra datasets have very different mean UMIs/genes, so be careful with intreptations (make sure to plot datasets means or vln plots)

# Need to subset hypo.integrated for only those subclusters that are species-specific, then plot these things for only zeb or only ast, with lines
sp.sp <- paste(names(tables), sp.sp.subclusters, sep = "_")
Idents(hypo.integrated) <- "integrated_Subcluster"
zeb.sp.sp <- subset(hypo.integrated, idents = sp.sp[sp %in% "zebrafish"])
ast.sp.sp <- subset(hypo.integrated, idents = sp.sp[sp %in% "tetra"])

hypo.integrated.zeb[["percent.mt"]] <- PercentageFeatureSet(hypo.integrated.zeb, pattern = "^mt-")
zeb.sp.sp[["percent.mt"]] <- PercentageFeatureSet(zeb.sp.sp, pattern = "^mt-")

zeb.mt <- VlnPlot(zeb.sp.sp, features = c("percent.mt"), group.by = "integrated_Subcluster", ncol = 1, pt.size = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 8), axis.title.x = element_blank(), axis.title.y = element_text(size = 10))
zeb.mt <- zeb.mt + ylab("% mito") + geom_hline(yintercept = c(mean(hypo.integrated.zeb@meta.data$percent.mt), mean(hypo.integrated.zeb@meta.data$percent.mt) + sd(hypo.integrated.zeb@meta.data$percent.mt), mean(hypo.integrated.zeb@meta.data$percent.mt) - sd(hypo.integrated.zeb@meta.data$percent.mt)), linetype = c("solid", "11", "11"), size = 0.75)

zeb.feature <- VlnPlot(zeb.sp.sp, features = c("nFeature_RNA"), group.by = "integrated_Subcluster", ncol = 1, pt.size = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 8), axis.title.x = element_blank(), axis.title.y = element_text(size = 10))
zeb.feature <- zeb.feature + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.zeb@meta.data$nFeature_RNA), mean(hypo.integrated.zeb@meta.data$nFeature_RNA) + sd(hypo.integrated.zeb@meta.data$nFeature_RNA), mean(hypo.integrated.zeb@meta.data$nFeature_RNA) - sd(hypo.integrated.zeb@meta.data$nFeature_RNA)), linetype = c("solid", "11", "11"), size = 0.75)
ast.feature <- VlnPlot(ast.sp.sp, features = c("nFeature_RNA"), group.by = "integrated_Subcluster", ncol = 1, pt.size = 0) + theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 8), axis.title.x = element_blank(), axis.title.y = element_text(size = 10))
ast.feature <- ast.feature + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.ast@meta.data$nFeature_RNA), mean(hypo.integrated.ast@meta.data$nFeature_RNA) + sd(hypo.integrated.ast@meta.data$nFeature_RNA), mean(hypo.integrated.ast@meta.data$nFeature_RNA) - sd(hypo.integrated.ast@meta.data$nFeature_RNA)), linetype = c("solid", "11", "11"), size = 0.75)

zeb.count <- VlnPlot(zeb.sp.sp, features = c("nCount_RNA"), group.by = "integrated_Subcluster", ncol = 1, pt.size = 0) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_blank(), axis.title.y = element_text(size = 10))
zeb.count <- zeb.count + ylab("# of counts") + geom_hline(yintercept = c(mean(hypo.integrated.zeb@meta.data$nCount_RNA), mean(hypo.integrated.zeb@meta.data$nCount_RNA) + sd(hypo.integrated.zeb@meta.data$nCount_RNA), mean(hypo.integrated.zeb@meta.data$nCount_RNA) - sd(hypo.integrated.zeb@meta.data$nCount_RNA)), linetype = c("solid", "11", "11"), size = 0.75)
ast.count <- VlnPlot(ast.sp.sp, features = c("nCount_RNA"), group.by = "integrated_Subcluster", ncol = 1, pt.size = 0) + theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8), axis.title.x = element_blank(), axis.title.y = element_text(size = 10))
ast.count <- ast.count + ylab("# of counts") + geom_hline(yintercept = c(mean(hypo.integrated.ast@meta.data$nCount_RNA), mean(hypo.integrated.ast@meta.data$nCount_RNA) + sd(hypo.integrated.ast@meta.data$nCount_RNA), mean(hypo.integrated.ast@meta.data$nCount_RNA) - sd(hypo.integrated.ast@meta.data$nCount_RNA)), linetype = c("solid", "11", "11"), size = 0.75)

layout <- "
AA###
BBDDD
CCEEE"

dev.new()
zeb.mt + zeb.feature + zeb.count + ast.feature + ast.count + plot_layout(design = layout, guides = "collect", heights = unit(c(30), c('mm')), widths = unit(c(20,20,20,30,30), c('mm')))



# ggplot(zeb.sp.sp@meta.data, aes(y = percent.mt, x = integrated_Subcluster, colour = integrated_Subcluster)) + geom_jitter(size = 1) + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.zeb@meta.data$percent.mt), mean(hypo.integrated.zeb@meta.data$percent.mt) + sd(hypo.integrated.zeb@meta.data$percent.mt), mean(hypo.integrated.zeb@meta.data$percent.mt) - sd(hypo.integrated.zeb@meta.data$percent.mt)), linetype = c("solid", "dotted", "dotted"), size = 1.5)
# ggplot(zeb.sp.sp@meta.data, aes(y = nFeature_RNA, x = integrated_Subcluster, colour = integrated_Subcluster)) + geom_jitter(size = 1) + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.zeb@meta.data$nFeature_RNA), mean(hypo.integrated.zeb@meta.data$nFeature_RNA) + sd(hypo.integrated.zeb@meta.data$nFeature_RNA), mean(hypo.integrated.zeb@meta.data$nFeature_RNA) - sd(hypo.integrated.zeb@meta.data$nFeature_RNA)), linetype = c("solid", "dotted", "dotted"), size = 1.5)
# ggplot(zeb.sp.sp@meta.data, aes(y = nCount_RNA, x = integrated_Subcluster, colour = integrated_Subcluster)) + geom_jitter(size = 1) + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.zeb@meta.data$nFeature_RNA), mean(hypo.integrated.zeb@meta.data$nFeature_RNA) + sd(hypo.integrated.zeb@meta.data$nFeature_RNA), mean(hypo.integrated.zeb@meta.data$nFeature_RNA) - sd(hypo.integrated.zeb@meta.data$nFeature_RNA)), linetype = c("solid", "dotted", "dotted"), size = 1.5)
# 
# ggplot(ast.sp.sp@meta.data, aes(y = nFeature_RNA, x = integrated_Subcluster, colour = integrated_Subcluster)) + geom_jitter(size = 1) + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.ast@meta.data$nFeature_RNA), mean(hypo.integrated.ast@meta.data$nFeature_RNA) + sd(hypo.integrated.ast@meta.data$nFeature_RNA), mean(hypo.integrated.ast@meta.data$nFeature_RNA) - sd(hypo.integrated.ast@meta.data$nFeature_RNA)), linetype = c("solid", "dotted", "dotted"), size = 1.5)
# ggplot(ast.sp.sp@meta.data, aes(y = nCount_RNA, x = integrated_Subcluster, colour = integrated_Subcluster)) + geom_jitter(size = 1) + ylab("# of genes") + geom_hline(yintercept = c(mean(hypo.integrated.ast@meta.data$nFeature_RNA), mean(hypo.integrated.ast@meta.data$nFeature_RNA) + sd(hypo.integrated.ast@meta.data$nFeature_RNA), mean(hypo.integrated.ast@meta.data$nFeature_RNA) - sd(hypo.integrated.ast@meta.data$nFeature_RNA)), linetype = c("solid", "dotted", "dotted"), size = 1.5)

