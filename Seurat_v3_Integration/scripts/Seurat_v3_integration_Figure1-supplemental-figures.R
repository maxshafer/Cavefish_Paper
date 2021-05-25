library(Seurat)
library(dplyr)
library(ggplot2)
library(reshape2)
library(grid)
library(viridis)
library(patchwork)

# Load objects

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")

hypo.integrated <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/Hypo_integrated_127k_1500VFs_100Dims_vR.rds")
hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k_vR.rds")
hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k_vR.rds")


## Make Prop plots
# Zeb
prop.table <- table(hypo.zeb@meta.data$Cluster, hypo.zeb@meta.data$orig.ident)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$cell_type <- row.names(prop.table)
prop.table <- reshape::melt(prop.table)
prop.table$cell_type <- as.factor(prop.table$cell_type)
prop.plots.zeb <- ggplot(prop.table, aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity")
prop.plots.zeb <- prop.plots.zeb + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_blank()) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + ylab("Sample Cluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Origin (cave or surface)")) + scale_fill_manual(values = cols3)

# Ast
prop.table <- table(hypo.ast@meta.data$Cluster, hypo.ast@meta.data$orig.ident)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$cell_type <- row.names(prop.table)
prop.table <- reshape::melt(prop.table)
prop.table$cell_type <- as.factor(prop.table$cell_type)
prop.plots.ast <- ggplot(prop.table, aes(x=cell_type, y=value, fill=variable)) + geom_bar(stat="identity")
prop.plots.ast <- prop.plots.ast + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.x = element_text(size = 10), axis.title.y = element_blank(), plot.background = element_rect(fill = "transparent", color = NA), panel.background = element_blank()) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + ylab("Sample Cluster frequency") + coord_flip() + guides(fill = guide_legend(title = "Origin (cave or surface)"))  + scale_fill_manual(values = cols3)

# Plot them together
pdf("Figures/Zeb-Ast_prop-plots.pdf", width = 10, height = 10)
prop.plots.zeb + prop.plots.ast + plot_layout(nrow = 2, guides = "collect", height = unit(c(87), c("mm")), width = unit(c(28), c("mm")))
dev.off()

## Make Marker gene plots

# Zeb
gene.lists <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/marker_gene_lists.rds")
genes.to.plot.zeb <- lapply(gene.lists[[1]], function(x) row.names(x)[1:2])
zeb.markers <- DotPlot(hypo.zeb, features = rev(unique(unlist(genes.to.plot.zeb))), group.by = "Cluster", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank())

# Ast
gene.lists <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/marker_gene_lists.rds")
genes.to.plot.ast <- lapply(gene.lists[[1]], function(x) row.names(x)[1:2])
# Replace with another gene (#2 and 3 aren't the best, 5 is good)
genes.to.plot.ast$Neuronal_07 <- row.names(gene.lists[[1]]$Neuronal_07)[c(1,5)]
ast.markers <- DotPlot(hypo.ast, features = rev(unique(unlist(genes.to.plot.ast))), group.by = "Cluster", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5), axis.title = element_blank())

# Plot them together

pdf("Figures/Zeb-Ast_marker-genes_plots.pdf", width = 10, height = 10)
zeb.markers + ast.markers + plot_layout(ncol = 2, guides = "collect", height = unit(c(140), c("mm")), width = unit(c(63), c("mm")))
dev.off()

## Make NT Dotplots
# Zeb
Idents(hypo.zeb) <- "Cluster"
hypo.zeb.neuronal <- subset(hypo.zeb, ident = levels(Idents(hypo.zeb))[grep("Neuronal", levels(Idents(hypo.zeb)))])

nt.plot.zeb <- DotPlot(object = hypo.zeb.neuronal, features = rev(c("slc17a6b", "gad2")), group.by = "Subcluster", scale.max = 200) + theme(legend.position = "right") + RotatedAxis() + scale_color_viridis()
nt.plot.zeb <- nt.plot.zeb + theme(axis.text = element_text(size = 6), axis.title = element_blank())

# Ast
Idents(hypo.ast) <- "Cluster"
hypo.ast.neuronal <- subset(hypo.ast, ident = levels(Idents(hypo.ast))[grep("Neuronal", levels(Idents(hypo.ast)))])

nt.plot.ast <- DotPlot(object = hypo.ast.neuronal, features = rev(c("slc17a6a","gad1b")), group.by = "Subcluster", scale.max = 200) + theme(legend.position = "right") + RotatedAxis() + scale_color_viridis()
nt.plot.ast <- nt.plot.ast + theme(axis.text = element_text(size = 6), axis.title = element_blank())

# Integrated
Idents(hypo.integrated) <- "integrated_Cluster"
hypo.integrated.neuronal <- subset(hypo.integrated, ident = levels(Idents(hypo.integrated))[grep("Neuronal", levels(Idents(hypo.integrated)))])

nt.plot <- DotPlot(object = hypo.integrated.neuronal, features = rev(c("slc17a6a", "slc17a6b", "gad1b", "gad2")), group.by = "integrated_Subcluster", scale.max = 200) + theme(legend.position = "right") + RotatedAxis() + scale_color_viridis()
nt.plot <- nt.plot + theme(axis.text = element_text(size = 6), axis.title = element_blank())

# Plot them

pdf("Figures/Neurotransmitter_plots.pdf", width = 10, height = 10)
nt.plot.zeb + nt.plot.ast + nt.plot + plot_layout(ncol = 3, guides = "collect", height = unit(c(200), c("mm")), width = unit(c(20,20,25), c("mm"))) 
dev.off()

## Make Integrated figure
# Markers
gene.lists <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/drift_gene_lists_2.rds")
genes.to.plot <- lapply(gene.lists[[1]], function(x) row.names(x)[1:5])
marker.sub <- DotPlot(hypo.integrated, features = rev(unique(unlist(genes.to.plot))), group.by = "integrated_Cluster", scale.max = 200) + coord_flip() + scale_color_viridis() + theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 8, angle = 45, hjust = 1, vjust = 1), axis.title = element_blank())
marker.sub <- marker.sub + theme(axis.text = element_text(size = 6))

# Integrated prop plots
props <- table(hypo.integrated@meta.data$integrated_Cluster, hypo.integrated@meta.data$species)
type_total <- apply(props, 1, function(x) sum(x[1],x[2],x[3]))
type_total <- (type_total/sum(type_total)*1000)
props[,1] <- props[,1]/sum(props[,1])*100
props[,2] <- props[,2]/sum(props[,2])*100
props[,3] <- props[,3]/sum(props[,3])*100
props <- melt(props)
props$type_total <- rep(type_total, 3)

# Order by total population level (can change to another thing)
cols <- c("#FDE725FF", "#22A884FF", "#414487FF")
props$Var.1 <- factor(props$Var.1, levels = levels(hypo.integrated@meta.data$integrated_Cluster))
prop.plot <- ggplot(props, aes(x = value, y = Var.1, width = (type_total), fill = Var.2)) + geom_bar(stat = "identity", position = "fill", colour = "black", size = 0) + scale_fill_manual(values = cols) + theme(panel.spacing.y = unit(.25, "mm"), strip.text.y = element_blank(), strip.background = element_blank(), axis.ticks.length.y =  unit(.5, "cm"), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(size = 8), axis.title.x = element_text(size = 10), legend.title = element_blank()) + xlab("Prop. of total species cells") + ylab(element_blank()) + geom_vline(xintercept = c(0.3333, 0.6666), linetype = "dashed", color = "black") + facet_grid(rows = vars(Var.1), scales = "free", space = "free")
prop.plot <- prop.plot + theme(axis.title.x = element_blank())


# Put together
pdf("Figures/Integration_markers_prop_plots.pdf", width = 10, height = 10)
marker.sub + prop.plot + plot_layout(ncol = 2, guides = "collect", height = unit(c(180), c("mm")), width = unit(c(80,20), c("mm")))
dev.off()