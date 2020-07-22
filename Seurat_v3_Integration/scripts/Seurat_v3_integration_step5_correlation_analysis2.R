library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(Matrix)
library(dplyr)
library(data.table)
library(purrr)
library(reshape2)
library(datapasta)
library(corrgram)
library(corrplot)
library(RColorBrewer)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

# Load object and marker genes
load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

### Load averaged data

norm.cluster <- readRDS("Normed_expression_data.rds")

str(norm.cluster, max.level = 2)

List of 5
 $ Type2                    :List of 4
  ..$ hypo.zeb    :List of 15
  ..$ hypo.ast    :List of 15
  ..$ hypo.surface:List of 15
  ..$ hypo.cave   :List of 15
 $ Subtype                  :List of 4
  ..$ hypo.zeb    :List of 41
  ..$ hypo.ast    :List of 36
  ..$ hypo.surface:List of 36
  ..$ hypo.cave   :List of 36
 $ SubclusterType           :List of 4
  ..$ hypo.zeb    :List of 135
  ..$ hypo.ast    :List of 171
  ..$ hypo.surface:List of 169
  ..$ hypo.cave   :List of 170
 $ integrated_Subtype       :List of 4
  ..$ hypo.zeb    :List of 26
  ..$ hypo.ast    :List of 26
  ..$ hypo.surface:List of 26
  ..$ hypo.cave   :List of 26
 $ integrated_SubclusterType:List of 4
  ..$ hypo.zeb    :List of 181
  ..$ hypo.ast    :List of 189
  ..$ hypo.surface:List of 186
  ..$ hypo.cave   :List of 187
  
  
  
  
## Create matrices for hypo.zeb and hypo.ast for all integrated Subtypes and SubclusterTypes

zeb.avg <- Reduce(cbind, norm.cluster[[5]][[1]])
colnames(zeb.avg) <- names(norm.cluster[[5]][[1]])

ast.avg <- Reduce(cbind, norm.cluster[[5]][[2]])
colnames(ast.avg) <- names(norm.cluster[[5]][[2]])

surface.avg <- Reduce(cbind, norm.cluster[[5]][[3]])
colnames(surface.avg) <- names(norm.cluster[[5]][[3]])

cave.avg <- Reduce(cbind, norm.cluster[[5]][[4]])
colnames(cave.avg) <- names(norm.cluster[[5]][[4]])


intersect.genes <- intersect(row.names(zeb.avg), row.names(ast.avg))

zeb.avg <- zeb.avg[intersect.genes,]
ast.avg <- ast.avg[intersect.genes,]
surface.avg <- surface.avg[intersect.genes,]
cave.avg <- cave.avg[intersect.genes,]

zeb.ast <- cor(zeb.avg, ast.avg) # row x column
zeb.surface <- cor(zeb.avg, surface.avg)
zeb.cave <- cor(zeb.avg, cave.avg)

zeb.diff <- zeb.surface[,intersect(colnames(zeb.surface), colnames(zeb.cave))] - zeb.cave[,intersect(colnames(zeb.surface), colnames(zeb.cave))]

## Add new Type annotation for integrated cell types

integrated_Subtypes <- levels(hypo.integrated@meta.data$integrated_Subtype)
Type2 <- c("GABAergic", "Endothelial", "Ependymal", "Glutamatergic", "Blood-Immune", "Glutamatergic", "Blood-Immune", "GABAergic", "GABAergic", "Glutamatergic", "Glutamatergic", "Progenitors", "GABAergic", "Glutamatergic", "GABAergic", "Glutamatergic", "Blood-Immune", "GABAergic", "GABAergic", "Glutamatergic", "Oligodendrocytes", "OPCs", "Lymphatic", "Ciliated", "Blood-Immune")
Type3 <- c("Neuronal", "Endothelial", "Glial", "Neuronal", "Blood-Immune", "Neuronal", "Blood-Immune", "Neuronal", "Neuronal", "Neuronal", "Neuronal", "Glial", "Neuronal", "Neuronal", "Neuronal", "Neuronal", "Blood-Immune", "Neuronal", "Neuronal", "Neuronal", "Glial", "Glial", "Lymphatic", "Glial", "Blood-Immune")

hypo.integrated@meta.data$Type2 <- as.character(plyr::mapvalues(hypo.integrated@meta.data$integrated_Subtype, integrated_Subtypes, Type2))
hypo.integrated@meta.data$Type3 <- as.character(plyr::mapvalues(hypo.integrated@meta.data$integrated_Subtype, integrated_Subtypes, Type3))

RowSideColors <- hypo.integrated@meta.data$Type3[match(rownames(zeb.ast), hypo.integrated@meta.data$integrated_SubclusterType)] # Should be based on zeb annotations
ColSideColors <- hypo.integrated@meta.data$Type3[match(colnames(zeb.ast), hypo.integrated@meta.data$integrated_SubclusterType)]

RowSideColors <- factor(RowSideColors, levels = unique(RowSideColors))
ColSideColors <- factor(ColSideColors, levels = unique(RowSideColors))

colours <- hue_pal()(length(levels(RowSideColors)))
names(colours) <- levels(RowSideColors)

RowSideColors <- as.character(plyr::mapvalues(RowSideColors, levels(RowSideColors), colours))
ColSideColors <- as.character(plyr::mapvalues(ColSideColors, levels(ColSideColors), colours))


heatmap.2(zeb.ast, trace = "none", Rowv = F, Colv = F, RowSideColors = RowSideColors, ColSideColors = ColSideColors)
heatmap.2(zeb.ast, trace = "none", RowSideColors = RowSideColors, ColSideColors = ColSideColors)
par(lend = 1)
legend(x = .88, y = 1.03, legend = names(colours), col = colours, lty = 1, lwd = 10)








zeb.ast <- melt(zeb.ast)
zeb.surface <- melt(zeb.surface)
zeb.cave <- melt(zeb.cave)
zeb.diff <- melt(zeb.diff)

zeb.ast$subtype.zeb <- hypo.integrated.zeb@meta.data$integrated_Subtype[match(zeb.ast$Var1, hypo.integrated.zeb@meta.data$integrated_SubclusterType)]
zeb.ast$subtype.ast <- hypo.integrated.ast@meta.data$integrated_Subtype[match(zeb.ast$Var2, hypo.integrated.ast@meta.data$integrated_SubclusterType)]

ggplot(zeb.ast, aes(x= Var1, y = Var2, fill = value)) + geom_tile(color = "transparent") + theme(axis.text = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_text(size = 10), panel.background = element_blank()) + scale_fill_viridis(option = "D") + ylab("Integrated subclusters - D. rerio") + xlab("Integrated subclusters - A. mexicanus")

ggplot(zeb.ast[intersect(grep("GABA_1", zeb.ast$Var1), grep("GABA_1", zeb.ast$Var2)),], aes(x= Var1, y = Var2, fill = value)) + geom_tile(color = "transparent") + theme(axis.text = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_text(size = 10), panel.background = element_blank()) + scale_fill_viridis(option = "D") + ylab("Integrated subclusters - D. rerio") + xlab("Integrated subclusters - A. mexicanus")

ggplot(zeb.ast[intersect(grep("GABA_[3-5]", zeb.ast$Var1), grep("GABA_[3-5]", zeb.ast$Var2)),], aes(x= Var1, y = Var2, fill = value)) + geom_tile(color = "transparent") + theme(axis.text = element_text(size = 6), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.title = element_text(size = 10), panel.background = element_blank()) + scale_fill_viridis(option = "D") + ylab("Integrated subclusters - D. rerio") + xlab("Integrated subclusters - A. mexicanus")

DotPlot(subset(hypo.integrated, idents = "GABA_1"), features = c("galn", "crabp1a", "otpb", "oxt", "sst1.1", "sst1.2", "pmchl", "creb3l1", "npvf", "hcrt", "lhx9"), group.by = "integrated_SubclusterType")
  
# Rbind the lists to use below
# Need to rbind each of the 2nd level lists

norm.cluster[[1]][[1]] <- do.call("cbind", norm.cluster[[1]][[1]])
colnames(norm.cluster[[1]][[1]]) <- 

  
# Get lists of peptides etc from zebrafish mine

markers.Subtype <- readRDS("Zeb_Ast_markers.Subtypes.filtered.rds")

markers.Subtype <- markers.Subtype[c(grep("p_val", markers.Subtype))]

for(i in 1:length(markers.Subtype)) {
	markers.Subtype[[i]]$gene <- rownames(markers.Subtype[[i]])
	markers.Subtype[[i]]$cluster <- names(markers.Subtype)[[i]]
}

markers.Subtype <- do.call(rbind, markers.Subtype)


# Get Go lists and use only those that are marker genes

go_lists <- list.files()[grep("GO", list.files())][3:8]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
go_lists <- lapply(go_lists, function(x) x[x %in% unique(markers.Subtype$gene)])

names(go_lists) <- list.files()[grep("GO", list.files())][3:8]
names(go_lists)

test <- c(as.character(unique(go_lists[[2]])), as.character(unique(go_lists[[5]])))






# Can we look at correlation amoung these genes?

matrix.zeb <- t(as.matrix(norm.cluster.zeb[test,]))
matrix.cor.zeb <- cor(matrix.zeb)
matrix.cor.zeb[is.na(matrix.cor.zeb)] <- 0

hr <- hclust(as.dist(1-matrix.cor.zeb), method = "complete")
order <- hr[[4]][hr[[3]]]
order.tf <- order[order %in% c(as.character(unique(go_lists[[2]])))]
order.np <- order[order %in% c(as.character(unique(go_lists[[5]])))]

matrix.cor.zeb.2 <- matrix.cor.zeb[order.np,order.tf]

corrplot(matrix.cor.zeb.2, type = "full", method = "color", order = "original", tl.col = "#666666", tl.srt = 45, tl.cex = 0.5, diag = T, addgrid.col = NA, cl.ratio = 0.07)

# Dot plots of top n correlated genes

gene <- "rims1a"

DotPlot(object = hypo.zeb, genes.plot = unique(c(names(tail(sort(matrix.cor.zeb.2[gene,]), n = 5)), gene)), plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType", do.return = T) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)) + coord_flip()







# Can we look at correlation amoung these genes?

matrix.ast <- t(as.matrix(norm.cluster.ast[test,]))
matrix.cor.ast <- cor(matrix.ast)
matrix.cor.ast[is.na(matrix.cor.ast)] <- 0

hr <- hclust(as.dist(1-matrix.cor.ast), method = "complete")
order <- hr[[4]][hr[[3]]]
order.tf <- order[order %in% c(as.character(unique(go_lists[[2]])))]
order.np <- order[order %in% c(as.character(unique(go_lists[[5]])))]

matrix.cor.ast.2 <- matrix.cor.ast[order.np,order.tf]

corrplot(matrix.cor.ast.2, type = "full", method = "color", order = "original", tl.col = "#666666", tl.srt = 45, tl.cex = 0.5, diag = T, addgrid.col = NA, cl.ratio = 0.07)

# Dot plots of top n correlated genes

gene <- "si:dkey-175g6.2"

DotPlot(object = hypo.ast, genes.plot = unique(c(names(tail(sort(matrix.cor.ast.2[gene,]), n = 5)), gene)), plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType", do.return = T) + theme(legend.position = "none", axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8)) + coord_flip()







hr <- hclust(as.dist(1-matrix.cor.zeb), method = "complete")

hc <- hclust(as.dist(1-t(matrix.cor.zeb)), method = "complete")

order <- hr[[4]][hr[[3]]]



matrix.cor.3 <- matrix.cor[order, order]
matrix.cor.3 <- matrix.cor.3[as.character(unique(go_lists[[3]])),]
matrix.cor.3 <- matrix.cor.3[,as.character(unique(go_lists[[4]]))]




DotPlot(object = hypo.ast, genes.plot = c(rownames(hypo.ast@data)[grep("pou", rownames(hypo.ast@data))], "hcrt"), plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType") + theme(legend.position = "right") + coord_flip()



DotPlot(object = hypo, genes.plot = rev(c("ltb4r", "mybl2b", "ets2", "nr4a2b", "bcl10", "nfkb2", "cebpa", "mafba", "rxrga", "atf4a", "nr4a3", "irak3", "spi1a", "bach1b", "irf8", "irf5", "spi1b", "batf3", "ikzf1", "rel", "cebpb", "irf7")), plot.legend = TRUE, x.lab.rot = TRUE, group.by = "SubclusterType") + theme(legend.position = "right") + coord_flip()



