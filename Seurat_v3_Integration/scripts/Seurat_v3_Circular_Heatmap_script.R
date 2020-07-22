library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggdendro)
library(phylogram)
library(ggtree)
library(grid)
library(ggmap)
library(dendextend)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

# load objects

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")
hypo <- hypo.integrated
cols3 <- c("#FDE725FF", "#22A884FF", "#414487FF")

# Plot the dendrogram

dendrograms <- readRDS("Zeb_Ast_dendrograms_integrated.rds")
#dendrograms <- lapply(dendrograms, function(x) as.dendrogram(x))

offsets <- rev(c(1, 9, 7, 7, 6, 3, 7, 5, 3, 8, 3, 5, 8, 4, 8, 6, 7, 4, 2, 3, 2, 3, 13, 5, 4, 5, 7, 6, 3, 1, 5, 1, 6, 1, 2, 3)*2)
offsets <- offsets + 10

dendro.plot <- ggtree(as.phylo(as.dendrogram(dendrograms[[1]])), branch.length = "none", layout = "circular") + theme_tree(bgcolor = NA)
# dendro.plot <- rotate_tree(dendro.plot, 80) + theme(plot.background = element_blank())
# dendro.plot <- dendro.plot + geom_tiplab2(aes(angle=angle),offset = offsets)

# orders <- unlist(lapply(dendrograms[[1]], function(x) labels(x)))

# Proportion figures for sex and orig.ident
# Make tables of cell type proportions
# Need to make matrix, for Subtype x Subcluster # (and fill in blanks, in reverse)

# 3 colours
prop.table <- table(hypo@meta.data$integrated_SubclusterType, hypo@meta.data$species)
prop.table <- prop.table/colSums(prop.table)

prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))

prop.table$cell_type <- row.names(prop.table)
index <- unique(hypo@meta.data[, c("integrated_SubclusterType", "integrated_Subtype")])
prop.table$Subtype <- index[match(prop.table$cell_type, index[,1]), 2]

hypo@meta.data$res.subclustertype <- unlist(lapply(strsplit(as.character(hypo@meta.data$integrated_SubclusterType), "_"), function(x) unlist(x)[length(unlist(x))]))

index <- unique(hypo@meta.data[, c("integrated_SubclusterType", "res.subclustertype")])
prop.table$res.0.4 <- index[match(prop.table$cell_type, index[,1]), 2]

# If 3 colours
prop.table$colour <- rgb(1-prop.table$astyanax_cave, 1-prop.table$astyanax_surface, 0, 1)

# prop.table <- prop.table[,c(5,1,2,3,4)]
prop.table <- melt(prop.table)
prop.table <- prop.table[order(match(prop.table$Subtype, orders)),]

prop.table$cell_type <- factor(prop.table$cell_type, levels = rev(unique(prop.table$cell_type)))

## Do it myself, without the ggtree wrapper - MUCH BETTER

df <- dendro.plot$data
width = 3

prop.table.filt <- prop.table[prop.table$variable == "zebrafish",]
prop.table.filt$y <- df$y[match(prop.table.filt$Subtype, df$label)]
prop.table.filt$col <- ifelse(prop.table.filt$value < 0.1 | prop.table.filt$value > 0.9, 1, 0)
prop.table.filt$size <- ifelse(prop.table.filt$value < 0.1 | prop.table.filt$value > 0.9, 0, 1)
prop.table.filt$x <- max(df$x) + 3 + as.numeric(prop.table.filt$res.0.4) * width

# cols <- c("#FDE725FF", "#22A884FF")

# If you want to do 3 colour graph
p.3 <- dendro.plot + geom_tile(data = prop.table.filt, aes(x, y), width = width, inherit.aes = FALSE, fill = prop.table.filt$colour) + geom_tile(data = prop.table.filt, aes(x, y,color = col, size = size), fill = NA, width = width, inherit.aes = FALSE) + scale_color_gradient(low = NA, high = "black", na.value = "white")

# If you want to do a 2 colour graph
p.2 <- dendro.plot + geom_tile(data = prop.table.filt, aes(x, y, fill = value), width = width, inherit.aes = FALSE)  + scale_fill_viridis(option = "D", direction = -1) + geom_tile(data = prop.table.filt, aes(x, y,color = col, size = size), fill = NA, width = width, inherit.aes = FALSE) + scale_color_gradient(low = NA, high = "red", na.value = "white")

# + scale_fill_gradient(low = "#ef8a62", high = "#67a9cf", na.value = NA)
# + scale_fill_gradient(low = "#FDE725FF", high = "#414487FF", na.value = NA)

## Add text for cell type labels and for cluster numbers

mapping <- data.frame(label = prop.table.filt$res.0.4, x = prop.table.filt$x)
mapping$y <- prop.table.filt$y

mapping2 <- unique(data.frame(label = prop.table.filt$Subtype, y = prop.table.filt$y))
mapping2$x <- apply(mapping2, 1, function(x) max(prop.table.filt$x[prop.table.filt$Subtype == x[1]]) + 5)
mapping2$hjust <- c(rep(1, 13), rep(0, 12))
mapping2$angle <- apply(mapping2, 1, function(x) max(dendro.plot$data$angle[dendro.plot$data$label == x[1]], na.rm = T))
mapping3 <- mapping2[order(mapping2$angle),]
mapping3$angle <- sapply(seq_along(c(1:nrow(mapping3))), function(x) min(mapping3$angle) + (x-1)*(max(mapping3$angle) - min(mapping3$angle))/nrow(mapping3))
mapping2 <- mapping3[row.names(mapping2),]
mapping2$angle[1:13] <- mapping2$angle[1:13] + 180
mapping2$angle <- mapping2$angle + 83

new.row <- data.frame(label = "Subclusters", y = 26, x = mean(as.numeric(mapping2$x)), angle = 80, hjust = .75)

p2 <- p.2 + theme(legend.position = "right", legend.title = element_blank()) + guides(size = FALSE, color = FALSE) + geom_text(data = mapping2, aes(x = x, y = y, label = label), angle = mapping2$angle, color = "black", size = 6/2.856, inherit.aes = FALSE, hjust = mapping2$hjust) + geom_text(data = mapping, aes(x = x, y = y, label = label), color = "white", size = 6/2.856, inherit.aes = FALSE, hjust = 0.5) + geom_text(data = new.row, aes(x = x, y = y, label = label), angle = new.row$angle, color = "black", size = 10/2.856, inherit.aes = FALSE, hjust = new.row$hjust, vjust = .8)

p3 <- p.3 + theme(legend.position = "right", legend.title = element_blank()) + guides(size = FALSE, color = FALSE) + geom_text(data = mapping2, aes(x = x, y = y, label = label), angle = mapping2$angle, color = "black", size = 6/2.856, inherit.aes = FALSE, hjust = mapping2$hjust) + geom_text(data = mapping, aes(x = x, y = y, label = label), color = "white", size = 6/2.856, inherit.aes = FALSE, hjust = 0.5) + geom_text(data = new.row, aes(x = x, y = y, label = label), angle = new.row$angle, color = "black", size = 10/2.856, inherit.aes = FALSE, hjust = new.row$hjust, vjust = .8)


p2 <- rotate_tree(p2, 80) + theme(plot.background = element_blank())

p3 <- rotate_tree(p3, 80) + theme(plot.background = element_blank())




## Make legend for p3 plot

line <- 6
col <- 6
red <- sort(rep(c(0,0.2,0.4,0.6,0.8,1),col))
green <- rep(c(0,0.2,0.4,0.6,0.8,1),line)
num <- 0


num <- num+1
colors <- rgb(red,green,i,1)
colors <- colors[c(1:6,7:11,13:16,19:21,25:26,31)]
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), axes = FALSE, xlab = "", ylab = "" )
mtext("% Zebrafish vs Astyanax surface", side=3 , line=0.15 , col="black" , font=2)
mtext("% Zebrafish vs Astyanax cave", side=2 , line=0.15 , col="black" , font=2)

rect(xleft, ybottom, xright, ytop, border = "light gray", col = colors)

xleft <- rep((0:(col - 1)/col),line)
xleft <- xleft[c(1:6,7:11,13:16,19:21,25:26,31)]

ybottom <- sort(rep((0:(line - 1)/line),col),decreasing=T)
ybottom <- ybottom[c(1:6,7:11,13:16,19:21,25:26,31)]

xright <- rep((1:col/col),line)
xleft <- xright[c(1:6,7:11,13:16,19:21,25:26,31)]

ytop <- sort(rep((1:line/line),col),decreasing=T)
ytop <- ytop[c(1:6,7:11,13:16,19:21,25:26,31)]

rect( rep((0:(col - 1)/col),line) ,  sort(rep((0:(line - 1)/line),col),decreasing=T), rep((1:col/col),line), sort(rep((1:line/line),col),decreasing=T),  border = "light gray" , col=colors)


axis(2 , at=c(17,14,11,8,5,2)/18-0.035 , labels=c("0","0.2","0.4","0.6","0.8","1") , tick=F , lty=6 , pos=0.01)
axis(3 , at=c(1.5 , 3.5 , 5.5 , 7.5 , 9.5 , 11.5)/12-0.045 , labels=c("0","0.2","0.4","0.6","0.8","1") , tick=F , pos=-0.15)

