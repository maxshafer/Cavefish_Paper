library(alluvial)
library(Seurat, lib.loc="/Users/maxwellshafer/Library/R/3.5/Dev")
library(reshape)
library(scales)
library(ggalluvial)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

load("Hypo_integrated_130k_1500VFs_100Dims_v2.Robj")


## Using alluvial package

table <- table(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "SubclusterType"], hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "integrated_SubclusterType"])
table <- melt(table)
col <- c(hue_pal()(length(unique(table$Var.2))))
names(table) <- c("Astyanax_SubclusterType", "integrated_SubclusterType", "value")

png("Figures/river_test_ast.png", height = 50, width = 10, res = 500, units = "in")
alluvial(table[,1:2], freq = table$value, col = col, border = col, hide = table$value < 100)
dev.off()

table <- table(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish", "integrated_SubclusterType"], hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "zebrafish", "SubclusterType"])
table <- melt(table)
col <- c(hue_pal()(length(unique(table$Var.1))))
names(table) <- c("integrated_SubclusterType", "Danio_SubclusterType", "value")

png("Figures/river_test_zeb.png", height = 50, width = 10, res = 500, units = "in")
alluvial(table[,1:2], freq = table$value, col = col, border = col, hide = table$value < 100)
dev.off()


## ggplot method, slower then alluvial

hypo.integrated@meta.data$zeb.Subtype <- hypo.integrated@meta.data$Subtype
hypo.integrated@meta.data$zeb.Subtype[hypo.integrated@meta.data$species.2 == "astyanax"] <- "astyanax"
hypo.integrated@meta.data$ast.Subtype <- hypo.integrated@meta.data$Subtype
hypo.integrated@meta.data$ast.Subtype[hypo.integrated@meta.data$species.2 == "zebrafish"] <- "zebrafish"

table <- table(hypo.integrated@meta.data[, "zeb.Subtype"], hypo.integrated@meta.data[, "integrated_Subtype"], hypo.integrated@meta.data[, "ast.Subtype"])

table <- melt(table)
col <- c(hue_pal()(length(unique(table$Var.1))))

table$value <- ifelse(table$value < 1000, 0, table$value)

ggplot(as.data.frame(table), aes(y= value, axis1 = Var.1, axis2 = Var.2, axis3 = Var.3)) + geom_alluvium(aes(fill = Var.1), width = 1/12) + geom_stratum(width = 1/12, fill = "black", color = "grey") + geom_label(stat = "stratum", label.strata = T) + scale_x_discrete(limits = c("integrated_Subtype", "Astyanax_Subtype"), expand = c(0.5, 0.5)) + guides(fill = F) + coord_flip()



table <- table(hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "integrated_Subtype"], hypo.integrated@meta.data[hypo.integrated@meta.data$species.2 == "astyanax", "Subtype"])
table <- melt(table)
col <- c(hue_pal()(length(unique(table$Var.1))))

table$value <- ifelse(table$value < 200, 0, table$value)

ggplot(as.data.frame(table), aes(y= value, axis1 = Var.1, axis2 = Var.2)) + geom_alluvium(aes(fill = Var.1), width = 1/12) + geom_stratum(width = 1/12, fill = "black", color = "grey") + geom_label(stat = "stratum", label.strata = T) + scale_x_discrete(limits = c("integrated_Subtype", "Astyanax_Subtype"), expand = c(0.5, 0.5)) + guides(fill = F) + coord_flip()
+ scale_fill_brewer(type = "qual", palette = "Set1")





## Method to maybe have continuous diagram (plots every cell, so super slow...)

data <- hypo.integrated@meta.data[,c("species.2", "integrated_Subtype", "Subtype")]
data$cell <- row.names(hypo.integrated@meta.data)

data$Ast_Subtype <- ifelse(data$species.2 == "astyanax", paste(data$species.2, data$Subtype, sep = "_"), NA) 
data$Zeb_Subtype <- ifelse(data$species.2 == "zebrafish", paste(data$species.2, data$Subtype, sep = "_"), NA) 

data <- data[,c(2,4,5,6)]

test <- melt(data, id.vars = c("cell"))

p <- ggplot(test2, aes(x = variable, stratum = value, alluvium = cell, fill = value, label = value)) + geom_flow(stat = "alluvium", lode.guidance = "rightleft", color = "darkgrey") + geom_stratum()

png("Figures/river_test_2.png", height = 50, width = 30, res = 500, units = "in")
p
dev.off()







## Make a proportional bar plot (thickness prop to cluster size) - does not work so well

Idents(hypo.integrated) <- "species.2"

hypo <- SubsetData(hypo.integrated, ident.use = "astyanax")

# Make tables of cell type proportions
prop.table <- table(hypo@meta.data$Subtype, hypo@meta.data$species)

total <- apply(prop.table, 1, function(y) sum(y))
prop.table <- as.data.frame(t(apply(prop.table, 1, function(y) {y/sum(y)})))
prop.table$total <- total
prop.table$total <- prop.table$total/sum(prop.table$total)

prop.table$cell_type <- row.names(prop.table)
prop.table <- prop.table[,c(4,1,2,3)]
prop.table <- prop.table %>% gather(species_morph, freq, astyanax_cave:astyanax_surface)
prop.table$cell_type <- as.factor(prop.table$cell_type)

ggplot(prop.table, aes(x=cell_type, y=freq, fill=species_morph, width = total)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + theme(axis.line.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y = element_text(size = 10), axis.title.x = element_blank(), plot.background = element_rect(fill = "transparent", color = NA)) + geom_hline(yintercept = 0.5, color = "grey36", size = 1, linetype = "dashed") + scale_fill_manual(values = c("lightgoldenrod1", "springgreen4","skyblue2")) + ylab("Species morph \n cluster frequency")
