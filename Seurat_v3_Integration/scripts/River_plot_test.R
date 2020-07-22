library(Seurat)
library(dplyr)
library(riverplot)
library(RColorBrewer)
library(reshape)


Idents(hypo.integrated) <- "species.2"

hypo.int.ast <- SubsetData(hypo.integrated, ident.use= "astyanax")
hypo.int.dan <- SubsetData(hypo.integrated, ident.use = "zebrafish")

table <- table(hypo.int.ast@meta.data$integrated_SubclusterType, hypo.int.ast@meta.data$SubclusterType)

table <- melt(table)

names1 <- levels(table[,"Var.1"])
names2 <- levels(table[,"Var.2"])

var1 <- as.numeric(table[,"Var.1"])
var2 <- as.numeric(table[,"Var.2"])
var2 <- var2 + max(var1)

edges <- table
edges$Var.1 <- as.numeric(edges$Var.1)
edges$Var.2 <- as.numeric(edges$Var.2) + max(var1)

colnames(edges) <- c("N1", "N2", "Value")


x <- c(rep(1, length(unique(table$Var.1))), rep(2, length(unique(table$Var.2))))

labels <- c(as.character(unique(table$Var.1)), as.character(unique(table$Var.2)))

ID <- unique(c(var1, var2))

col <- c(hue_pal()(length(unique(table$Var.1))), hue_pal()(length(unique(table$Var.2))))

srt <- rep(0, length(col))

nodes <- data.frame(ID = ID, x = x, labels = labels, col = col, srt = srt)

river <- makeRiver(nodes, edges)

png("Figures/rivertest.png", width = 20, height = 100, units = "in", res = 250)
plot(river)
dev.off()



ids <- c("Subtype", "SubclusterType")

Idents(hypo.int.ast) <- "Subtype"
colours <- hue_pal()(length(levels(Idents(hypo.int.ast))))
names(colours) <- levels(Idents(hypo.int.ast))


my_colour_palette <- list()
for (i in 1:length(ids)) {
	hypo <- SetAllIdent(hypo, id = ids[[i]])
	colours <- hue_pal()(length(levels(hypo@ident)))
	names(colours) <- levels(hypo@ident)
	colours <- colours[dendrograms[[i]] %>% labels]
	my_colour_palette[[i]] <- colours
}





