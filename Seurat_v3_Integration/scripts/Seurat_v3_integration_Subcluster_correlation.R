library(Seurat)
library(Matrix)
library(dplyr)
library(scales)
library(dplyr)
library(tidyr)
library(reshape2)
library(magick)
library(ggimage)
library(png)
library(grid)
library(phytools)
library(ggtree)
library(ggplotify)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Seurat_v3_Integration/")

# load objects

load("/Volumes/BZ/Home/gizevo30/R_Projects/Shafer_Hypo/Shafer_Hypo_66k.Robj")
load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_64k.Robj")

normed.expression <- readRDS("Normed_expression_data.rds")

str(normed.expression, max.level = 2)

# ## Need to remove GABA_5, and rename GABA_6 in normed.expression[[4]] and [[5]]

# for(i in 1:length(normed.expression[[4]])) {
	# normed.expression[[4]][[i]] <- normed.expression[[4]][[i]][!grepl("GABA_5", names(normed.expression[[4]][[i]]))]
	# names(normed.expression[[4]][[i]])[grep("GABA_6", names(normed.expression[[4]][[i]]))] <- "GABA_5"
# }

# for(i in c(1,3,4)) {
	# normed.expression[[5]][[i]] <- normed.expression[[5]][[i]][!grepl("GABA_5", names(normed.expression[[5]][[i]]))]
	# names(normed.expression[[5]][[i]])[grep("GABA_6", names(normed.expression[[5]][[i]]))] <- c("GABA_5_0", "GABA_5_1", "GABA_5_2")
# }

# saveRDS(normed.expression, file = "Normed_expression_data.rds")


# cor(t(norm.cluster.zeb)[,2], t(norm.cluster.surface)[,2])

correlateExp <- function(zeb = zeb, ast = ast) {
	correlation <- list()
	for(i in 1:length(zeb)){
		genes <- intersect(row.names(zeb[[i]]), row.names(ast[[i]]))
		correlation[[i]] <- unlist(lapply(seq_along(zeb[[i]]), function(x) cor(zeb[[i]][genes,x], ast[[i]][genes,x])))
	}
	names(correlation) <- names(zeb)
	return(unlist(correlation))
}

corr.surface <- correlateExp(zeb = normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]], ast = normed.expression[["integrated_SubclusterType"]][["hypo.surface"]])

corr.cave <- correlateExp(zeb = normed.expression[["integrated_SubclusterType"]][["hypo.zeb"]], ast = normed.expression[["integrated_SubclusterType"]][["hypo.cave"]])

test <- data.frame(corr.surface, corr.cave)
test$Subclusters <- row.names(test)
test$Subtypes <- hypo.integrated@meta.data$integrated_Subtype[match(test$Subclusters, hypo.integrated@meta.data$integrated_SubclusterType)]
test$diff <- ifelse(abs(test$corr.surface - test$corr.cave) > 0.6, "yes", "no")

test$Subtypes <- factor(test$Subtypes, levels = c("Endothelial", "Erythrocytes", "Ciliated", "Ependymal", "Progenitors", "Oligodendrocyte_Precursor_Cells", "Oligodendrocytes", "GABA_0", "GABA_1", "GABA_2", "GABA_3", "GABA_4", "GABA_5", "Prdx1_Positive", "Glut_0", "Glut_1", "Glut_2", "Glut_3", "Glut_4", "Glut_5", "Glut_6", "Lymphatic", "Leucocytes", "Macrophages", "Microglia"))

ggplot(test, aes(x = corr.surface, y = corr.cave, color = Subtypes)) + geom_point() + xlab("D. rerio vs A. mexicanus - Surface") + ylab("D. rerio vs A. mexicanus - Cave") + theme(axis.text = element_text(size = 8), text = element_text(size = 10)) + geom_text_repel(data = test[test$diff == "yes",], aes(label = Subclusters), size = 2)


ggplot(test, aes(x = Subtypes, y = corr.surface - corr.cave, color = Subtypes)) + geom_jitter() + guides(color = F) + theme(axis.title = element_text(size = 10), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), axis.text = element_text(size = 8)) + xlab("") + ylab("Similary D. rerio \n vs A. mexicanus Surface") + ylim(-1,1) + geom_hline(yintercept = 0, linetype = "dashed", color = "black") + geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") + geom_hline(yintercept = -0.5, linetype = "dashed", color = "red")
