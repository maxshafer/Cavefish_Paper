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

load("Shafer_Hypo_66k.Robj")
hypo <- hypo.zeb

markers.SubclusterType.zeb <- readRDS("Shafer_Hypo_markers.SubclusterType.rds")
markers.SubclusterType.ast <- readRDS(file = "/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/species.subcluster.markers.rds")

normed.sc.zeb <- readRDS("Shafer_Hypo_normed_SubclusterType_matrix.rds")
normed.sc.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_normed_SubclusterType_matrix.rds")

# Get lists of peptides etc from zebrafish mine

go_lists <- list.files()[grep("GO", list.files())]
go_lists <- lapply(go_lists, function(x) read.csv(x, head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))
go_lists <- lapply(go_lists, function(x) x[x %in% markers.SubclusterType.zeb$gene])
go_lists <- lapply(go_lists, function(x) x[x %in% markers.SubclusterType.ast$gene])

names(go_lists) <- list.files()[grep("GO", list.files())]
names(go_lists)

test <- c(as.character(unique(go_lists[[2]])), as.character(unique(go_lists[[5]])))

# Can we look at correlation amoung these genes?

matrix.zeb <- t(as.matrix(normed.sc.zeb[test,]))
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

matrix.ast <- t(as.matrix(normed.sc.ast[test,]))
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



