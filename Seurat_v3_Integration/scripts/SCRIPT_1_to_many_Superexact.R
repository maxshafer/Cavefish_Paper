library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(reshape)
library(scales)
library(ggalluvial)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(SuperExactTest)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")

Idents(hypo.integrated) <- "integrated_Subtype"

source("~/Documents/Seurat_objects/R_script_load_objects.R")


## Get species specific cell types

prop.table <- table(hypo.integrated@meta.data$integrated_SubclusterType, hypo.integrated@meta.data$species)
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$Zeb_specific <- ifelse(prop.table$zebrafish > .9, "yes", "no")
prop.table$Ast_specific <- ifelse(prop.table$zebrafish < .1, "yes", "no")

zeb.names <- row.names(prop.table[prop.table$Zeb_specific == "yes",])
ast.names <- row.names(prop.table[prop.table$Ast_specific == "yes",])

species.specific <- c(zeb.names, ast.names)

#### THIS WORKS
#### FOR SUBCLUSTERS


table1 <- table(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish", "integrated_SubclusterType"], hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish","SubclusterType"])
table2 <- table(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "integrated_SubclusterType"], hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax","SubclusterType"])

## Generate a confusion matrix!! Makes more sense
# Zeb percentage per integrated, percentage of astyanax
# For each zeb row, make a new matrix

## Need to add columns to zeb.percentage, and rows to ast.percentage
## Maybe best to make this the same way, then 

zeb.percentage <- (table1/rowSums(table1))*100 # integrated_SubclusterType x SubclusterType
ast.percentage <- (table2/rowSums(table2))*100 # integrated_SubclusterType x SubclusterType

all.ids <- unique(hypo.integrated@meta.data$integrated_SubclusterType)

zeb.percentage.all <- do.call(rbind,lapply(all.ids, function(x) {
	if(x %in% rownames(zeb.percentage)){
		zeb.percentage[x,]
	} else {
		rep(0, ncol(zeb.percentage))
	}
}))
rownames(zeb.percentage.all) <- all.ids

ast.percentage.all <- do.call(rbind,lapply(all.ids, function(x) {
	if(x %in% rownames(ast.percentage)){
		ast.percentage[x,]
	} else {
		rep(0, ncol(ast.percentage))
	}
}))
rownames(ast.percentage.all) <- all.ids

# zeb.percentage.all[zeb.percentage.all < 10] <- 0
# ast.percentage.all[ast.percentage.all < 10] <- 0

zeb.percentage.all <- zeb.percentage.all[!(row.names(zeb.percentage.all) %in% species.specific),]
ast.percentage.all <- ast.percentage.all[!(row.names(ast.percentage.all) %in% species.specific),]

## Make Pie chart!

zeb.counts <- rowSums(ifelse(zeb.percentage.all > 10, as.numeric(1), as.numeric(0)), na.rm = T)
ast.counts <- rowSums(ifelse(ast.percentage.all > 10, as.numeric(1), as.numeric(0)), na.rm = T)

combined.counts <- data.frame(zeb.counts = zeb.counts, ast.counts = ast.counts)
counts <- plyr::count(combined.counts)
counts$x <- apply(counts, 1, function(x) paste(ifelse(x[1] < 3, x[1], "many"), ":", ifelse(x[2] < 3, x[2], "many"), sep = ""))
counts <- plyr::count(counts, vars = "x")
counts$x <- factor(counts$x, levels = counts$x[c(1,4,7,5,6,8,9,2,3)])
counts <- counts[c(1,4,7,5,6,8,9,2,3),]
counts$colour <- rev(c(rep("lightblue", 2), rep("yellow", 4), rep("red", 2), "black"))
counts$freq2 <- 2
counts$freq3 <- 1

pie.chart <- ggplot(counts, aes(x = freq3, y = freq, fill = x)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + scale_fill_manual(values = c("black", brewer.pal(8, "RdYlBu")))

pie.chart + geom_bar(data = counts, aes(x = freq2, y = freq, color = colour), width = .01, stat = "identity", fill = "transparent") + coord_polar("y", start = 0) + xlim(c(.5,2.1)) + scale_color_manual(values = (unique(counts$colour))) + theme(axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + labs(fill = "D. rerio : A. mexicanus") + geom_text(aes(y = rev(cumsum(freq) - 0.5*freq), label = freq), size = 4)


### Make lists of 1:1, 1:2 etc integrated subclusters

zeb.counts <- rowSums(ifelse(zeb.percentage.all > 10, as.numeric(1), as.numeric(0)), na.rm = T)
ast.counts <- rowSums(ifelse(ast.percentage.all > 10, as.numeric(1), as.numeric(0)), na.rm = T)

combined.counts <- data.frame(zeb.counts = zeb.counts, ast.counts = ast.counts)

combined.counts$x <- apply(combined.counts, 1, function(x) paste(ifelse(x[1] < 3, x[1], "many"), ":", ifelse(x[2] < 3, x[2], "many"), sep = ""))
combined.counts$SubclusterType <- row.names(combined.counts)

orthologies <- lapply(unique(combined.counts$x), function(x) combined.counts[combined.counts$x == x, "SubclusterType"])
names(orthologies) <- unique(combined.counts$x)

### Make list of zebrafish and astyanax subclusters associated with each integrated subcluster

zeb.orthologies <- lapply(row.names(zeb.percentage.all), function(x) colnames(zeb.percentage.all)[zeb.percentage.all[x,] > 10])
names(zeb.orthologies) <- row.names(zeb.percentage.all)

ast.orthologies <- lapply(row.names(ast.percentage.all), function(x) colnames(ast.percentage.all)[ast.percentage.all[x,] > 10])
names(ast.orthologies) <- row.names(ast.percentage.all)


## Subset and find HVGs
detach("package:Seurat", unload = T)
library(Seurat)

hypo.zeb <- SetAllIdent(hypo.zeb, id = "integrated_SubclusterType")

# zeb.orthologies.subset <- zeb.orthologies[c(orthologies[[2]], orthologies[[5]])]
zeb.orthologies.subset <- zeb.orthologies[c(orthologies[[2]], orthologies[[5]], orthologies[[6]], orthologies[[7]], orthologies[[8]], orthologies[[9]])]

subset.hvg.zeb <- list()
for (i in 1:length(zeb.orthologies.subset)) {
	subset <- SubsetData(hypo.zeb, ident.use = names(zeb.orthologies.subset)[[i]])
	subset <- SetAllIdent(subset, id = "SubclusterType")
	subset <- SubsetData(subset, ident.use = zeb.orthologies.subset[[i]])
	subset <- FindVariableGenes(subset, x.low.cutoff = 0.5, y.cutoff = 2.5, display.progress = F, do.plot = F)
	subset.hvg.zeb[[i]] <- subset@var.genes
}
names(subset.hvg.zeb) <- names(zeb.orthologies.subset)

hypo.ast <- SetAllIdent(hypo.ast, id = "integrated_SubclusterType")

# ast.orthologies.subset <- ast.orthologies[c(orthologies[[1]], orthologies[[4]])]
ast.orthologies.subset <- ast.orthologies[c(orthologies[[1]], orthologies[[4]], orthologies[[6]], orthologies[[7]], orthologies[[8]], orthologies[[9]])]

subset.hvg.ast <- list()
for (i in 1:length(ast.orthologies.subset)) {
	subset <- SubsetData(hypo.ast, ident.use = names(ast.orthologies.subset)[[i]])
	subset <- SetAllIdent(subset, id = "SubclusterType")
	subset <- SubsetData(subset, ident.use = ast.orthologies.subset[[i]])
	subset <- FindVariableGenes(subset, x.low.cutoff = 0.5, y.cutoff = 2.5, display.progress = F, do.plot = F)
	subset.hvg.ast[[i]] <- subset@var.genes
}
names(subset.hvg.ast) <- names(ast.orthologies.subset)

detach("package:Seurat", unload = T)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")

# SuperExactTest

total <- nrow(hypo.ast@data)

# Can't possible do this for all possible combinations, so split into smaller groups of integrated cell types
subset <- subset.hvg.zeb[grep("Prog", names(subset.hvg.zeb))]

# For ast
subset <- subset.hvg.ast[grep("Glut", names(subset.hvg.ast))]

res=supertest(subset, n=total)

plot(res, Layout="landscape", degree = 2:10, sort.by="size", keep=FALSE, show.elements=TRUE, elements.cex=0.7, show.fold.enrichment = TRUE, elements.list=subset(summary(res)$Table,Observed.Overlap <= 12), show.expected.overlap=TRUE,expected.overlap.style="hatchedBox", color.expected.overlap='red')

test <- summary(res)$Table




# hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish", "Subtype"] <- hypo.zeb@meta.data$Subtype[match(rownames(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish",]), rownames(hypo.zeb@meta.data))]

# hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "Subtype"] <- hypo.ast@meta.data$Subtype[match(rownames(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax",]), rownames(hypo.ast@meta.data))]

# hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish", "SubclusterType"] <- hypo.zeb@meta.data$SubclusterType[match(rownames(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish",]), rownames(hypo.zeb@meta.data))]

# hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "SubclusterType"] <- hypo.ast@meta.data$SubclusterType[match(rownames(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax",]), rownames(hypo.ast@meta.data))]




# #### OLD SCRIPT, for calling zeb in ast, and ast in zeb (not relative to integrated)

# # Set up data

# table1 <- table(hypo@meta.data[hypo@meta.data$species.2 == "zebrafish", "Subtype"], hypo@meta.data[hypo@meta.data$species.2 == "zebrafish","integrated_Subtype"])
# table2 <- table(hypo@meta.data[hypo@meta.data$species.2 == "astyanax", "integrated_Subtype"], hypo@meta.data[hypo@meta.data$species.2 == "astyanax","Subtype"])

# ## Generate a confusion matrix!! Makes more sense
# # Zeb percentage per integrated, percentage of astyanax
# # For each zeb row, make a new matrix

# zeb.percentage <- (table1/rowSums(table1))*100
# ast.percentage <- t((t(table2)/rowSums(t(table2)))*100)

# zeb.to.ast <- apply(zeb.percentage, 1, function(x) colSums(ast.percentage*x, na.rm = T)/sum(colSums(ast.percentage*x, na.rm = T))*100)

# # For reverse
# ast.to.zeb <- apply(t(ast.percentage), 1, function(x) colSums(t(zeb.percentage)*x)/sum(colSums(t(zeb.percentage)*x))*100)

# ggplot(melt(zeb.to.ast*t(ast.to.zeb)), aes(x = Var.1, y = Var.2, fill = value)) + geom_tile() + scale_fill_viridis(discrete = F, option = "A", begin = 1, end = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 5), axis.title = element_blank())


# ## Make zeb.to.ast and ast.to.zeb matrices

# zeb.to.ast <- apply(zeb.percentage.all, 2, function(x) (colSums(x*ast.percentage.all, na.rm = T)/sum(colSums(x*ast.percentage.all, na.rm = T)))*100)
# is.na(zeb.to.ast) <- 0
# zeb.to.ast[zeb.to.ast == "NaN"] <- 0

# # For reverse
# ast.to.zeb <- apply(ast.percentage.all, 2, function(x) colSums(zeb.percentage.all*x, na.rm = T)/sum(colSums(zeb.percentage.all*x, na.rm = T))*100)
# is.na(ast.to.zeb) <- 0
# ast.to.zeb[ast.to.zeb == "NaN"] <- 0

# ## Plot matrices with ggplot!

# matrix.plots <- list()

# matrix.plots[[1]] <- ggplot(melt(ast.to.zeb), aes(x = X1, y = X2, fill = value)) + geom_tile() + scale_fill_viridis(discrete = F, option = "A", begin = 1, end = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 5))

# matrix.plots[[2]] <- ggplot(melt(t(zeb.to.ast)), aes(x = X1, y = X2, fill = value)) + geom_tile() + scale_fill_viridis(discrete = F, option = "A", begin = 1, end = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 5))

# product <- (ast.to.zeb*t(zeb.to.ast))/100
# matrix.plots[[3]] <- ggplot(melt(product), aes(x = X1, y = X2, fill = value)) + geom_tile() + scale_fill_viridis(discrete = F, option = "A", begin = 1, end = 0) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), axis.text = element_text(size = 5))

# matrix.legend <- get_legend(matrix.plots[[1]])
# matrix.plots2 <- lapply(matrix.plots, function(x) x + theme(legend.position = "none", axis.title = element_blank()))

# matrices <- plot_grid(plotlist = matrix.plots2, ncol = 3)
# matrices <- plot_grid(matrices, matrix.legend, rel_widths = c(1,.05))


# ## Make Orthologous cell type/subcluster calls, and plot pie charts!

# library(plyr)

# zeb.counts <- unlist(lapply(seq_along(colnames(zeb.to.ast)), function(x) length(zeb.to.ast[,x][zeb.to.ast[,x] > 10])))
# ast.counts <- unlist(lapply(seq_along(colnames(ast.to.zeb)), function(x) length(ast.to.zeb[,x][ast.to.zeb[,x] > 10])))
# product.counts <- unlist(lapply(seq_along(colnames(product)), function(x) length(product[,x][product[,x] > 10])))

# zeb.counts <- count(zeb.counts)
# ast.counts <- count(ast.counts)
# product.counts <- count(product.counts)

# zeb.counts$x <- c("1:none", "1:1", "1:2", "1:many", "1:many", "1:many", "1:many")
# ast.counts$x <- c("1:none", "1:1", "1:2", "1:many", "1:many", "1:many", "1:many")
# product.counts$x <- c("1:none", "1:1", "1:2", "1:3")

# product.counts <- count(product.counts, vars = "x")
# zeb.counts <- count(zeb.counts, vars = "x")
# ast.counts <- count(ast.counts, vars = "x")

# product.counts$species <- "D. rerio * A. mexicanus"
# zeb.counts$species <- "D. rerio"
# ast.counts$species <- "A. mexicanus"

# zeb.counts$freq <- zeb.counts$freq/sum(zeb.counts$freq)*100
# ast.counts$freq <- ast.counts$freq/sum(ast.counts$freq)*100
# product.counts$freq <- product.counts$freq/sum(product.counts$freq)*100

# # Make pie charts

# pies <- list()
# pies[[1]] <- ggplot(zeb.counts, aes(x = "", y = freq, fill = x)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + ggtitle("D. rerio -> A. mexicanus")
# pies[[2]] <- ggplot(ast.counts, aes(x = "", y = freq, fill = x)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + ggtitle("A. mexicanus -> D. rerio")
# pies[[3]] <- ggplot(product.counts, aes(x = "", y = freq, fill = x)) + geom_bar(width = 1, stat = "identity") + coord_polar("y", start = 0) + ggtitle("Mutual relationships")

# pies2 <- lapply(pies, function(x) x + scale_fill_brewer(palette = "RdGy") + theme(axis.text = element_blank(), axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank()) + geom_text(aes(y = freq, label = percent(freq/100)), size = 5))
# legend <- get_legend(pies2[[1]])

# pies2 <- lapply(pies2, function(x) x + theme(legend.position = "none"))
# pies3 <- plot_grid(plotlist = pies2, ncol = 3)
# pies3 <- plot_grid(pies3, legend, rel_widths = c(1,.1))
# pies3

# plot_grid(plot_grid(pies3, ncol = 2), matrices, nrow = 2, rel_heights = c(1,3))


