library(Seurat)
library(Matrix)
library(dplyr)
library(datapasta)
library(corrgram)
library(corrplot)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

# Generate list of markers for subclusters (gene.lists[[5]] for zeb, gene.lists[[6]] for ast)

markers.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_markers.Subtype.rds")
markers.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_markers.Subtype.rds")

# Find genes which are markers for multiple subclusters

markers.zeb.sub <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_markers.SubclusterType.rds")
markers.zeb.sub <- markers.zeb.sub %>% group_by(SubclusterType) %>% top_n(-50, value)
markers.zeb.sub <- tibble(gene = as.character(markers.zeb.sub$gene))
markers.zeb.sub <- markers.zeb.sub %>% group_by(gene) %>% tally(sort = T)
markers.zeb.sub <- markers.zeb.sub[markers.zeb.sub$n > 1,]

markers.zeb.sub <- markers.zeb.sub[!(markers.zeb.sub$gene %in% markers.zeb$gene[markers.zeb$p_val_adj > 0.001]),]


markers.ast.sub <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_Hypo_markers.SubclusterType.rds")
markers.ast.sub <- markers.ast.sub %>% group_by(SubclusterType) %>% top_n(-50, value)
markers.ast.sub <- tibble(gene = as.character(markers.ast.sub$gene))
markers.ast.sub <- markers.ast.sub %>% group_by(gene) %>% tally(sort = T)
markers.ast.sub <- markers.ast.sub[markers.ast.sub$n > 1,]

markers.ast.sub <- markers.ast.sub[!(markers.ast.sub $gene %in% markers.ast$gene[markers.ast$p_val_adj > 0.001]),]


zeb.genes <- ggplot(markers.zeb.sub[1:30,], aes(x = gene, y = n)) + geom_bar(stat = "identity") + coord_flip() + theme(axis.title = element_blank())
ast.genes <- ggplot(markers.ast.sub[1:30,], aes(x = gene, y = n)) + geom_bar(stat = "identity") + coord_flip() + theme(axis.title = element_blank())

overlap <- intersect(markers.ast.sub$gene, markers.zeb.sub$gene)
zeb <- markers.zeb.sub[markers.zeb.sub$gene %in% overlap,]
ast <- markers.ast.sub[markers.ast.sub$gene %in% overlap,]
test <- full_join(zeb, ast, by = "gene")

scatter <- ggplot(test, aes(x = n.x, y = n.y)) + geom_jitter() + xlab("Re-use in Zebrafish") + ylab("Re-use in Astyanax")

png("Figures/Marker_gene_reuse_scatter.png", width = 12, height = 4, units = "in", res = 250)
plot_grid(scatter, zeb.genes, ast.genes, ncol = 3, rel_widths = c(4,3,4))
dev.off()

# ggplot(melt(test[1:30,]), aes(x = gene, y = value, fill = variable)) + geom_bar(stat = "identity", position = "dodge") + coord_flip()


png("Figures/Marker_gene_reuse_zeb.png", width = 12, height = 3, units = "in", res = 250)
FeaturePlot(object = hypo.zeb, no.axes = T, reduction.use = "tsne", features.plot = c("penkb", "ccka", "adcyap1b", "fosb"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .25, do.return = F, nCol = 4)
dev.off()

png("Figures/Marker_gene_reuse_ast.png", width = 12, height = 3, units = "in", res = 250)
FeaturePlot(object = hypo.ast, no.axes = T, reduction.use = "FItSNE", features.plot = c("penkb", "ccka", "adcyap1b", "fosb"), cols.use = c("grey85", "blue"), min.cutoff = "q9", pt.size = .5, do.return = F, nCol = 4)
dev.off()


DotPlot(hypo.zeb, genes.plot = c("penkb", "ccka", "adcyap1b", "fosb"), group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE)

DotPlot(hypo.zeb, genes.plot = markers.zeb.sub$gene[1:30], group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE)
DotPlot(hypo.ast, genes.plot = markers.ast.sub$gene[1:30], group.by = "SubclusterType", plot.legend = TRUE, x.lab.rot = TRUE)



# # Can we look at correlation amoung these genes?

# markers.zeb.sub.cut <- markers.zeb.sub[markers.zeb.sub$n > 5,]

# matrix <- t(as.matrix(hypo.zeb@data[markers.zeb.sub.cut$gene,]))

# matrix.cor <- cor(matrix)

# corrplot(matrix.cor, type = "full", method = "color", order = "hclust", tl.col = "#666666", tl.srt = 45, tl.cex = 0.5, col = brewer.pal(n=15, name="PuOr"), diag = F, addgrid.col = NA, cl.ratio = 0.07)

