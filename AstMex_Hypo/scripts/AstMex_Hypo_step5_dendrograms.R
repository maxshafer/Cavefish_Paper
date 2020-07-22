library(Seurat)
library(tidyr)
library(dendextend)
library(dplyr)
library(ape)
library(circlize)
library(ggplot2)
library(scales)

# Load seurat object (astmex)

load("/Volumes/BZ/Home/gizevo30/R_Projects/AstMex_Hypo/AstMex_64k.Robj")

###################################################
# Make dendrograms
###################################################

# Extract dendrograms from hypo object (4 of them)

dendrograms <- hypo@cluster.tree
names(dendrograms)
# [1] "Subtype"            "SubclusterType"     "species_subtype"    "species_subcluster"

dendrograms <- lapply(dendrograms, function(x) as.dendrogram(x))

# Make colour palettes

ids <- names(dendrograms)

my_colour_palette <- list()
for (i in 1:length(ids)) {
	hypo <- SetAllIdent(hypo, id = ids[[i]])
	colours <- hue_pal()(length(levels(hypo@ident)))
	names(colours) <- levels(hypo@ident)
	colours <- colours[dendrograms[[i]] %>% labels]
	my_colour_palette[[i]] <- colours
}

# Prepare label colours (for species specific)

label.colours <- list()
label.colours <- lapply(dendrograms, function(x) labels(x))

label.colours[[1]] <- rep("black", length(label.colours[[1]]))
label.colours[[2]] <- rep("black", length(label.colours[[2]]))
label.colours[[3]][grep("astyanax_cave", label.colours[[3]])] = "lightgoldenrod1"
label.colours[[3]][grep("astyanax_surface", label.colours[[3]])] = "springgreen4"
label.colours[[4]][grep("astyanax_cave", label.colours[[4]])] = "lightgoldenrod1"
label.colours[[4]][grep("astyanax_surface", label.colours[[4]])] = "springgreen4"
label.colours[[5]][grep("choy_surface", label.colours[[5]])] = "springgreen4"
label.colours[[5]][grep("molino_cave", label.colours[[5]])] = "goldenrod1"
label.colours[[5]][grep("pachon_cave", label.colours[[5]])] = "lightgoldenrod1"
label.colours[[5]][grep("tinaja_cave", label.colours[[5]])] = "darkorange1"
label.colours[[6]][grep("choy_surface", label.colours[[6]])] = "springgreen4"
label.colours[[6]][grep("molino_cave", label.colours[[6]])] = "goldenrod1"
label.colours[[6]][grep("pachon_cave", label.colours[[6]])] = "lightgoldenrod1"
label.colours[[6]][grep("tinaja_cave", label.colours[[6]])] = "darkorange1"


# Prepare node labels

node.labels <- list()
node.labels <- lapply(dendrograms, function(x) labels(x))

node.labels[[1]] <- rep(19, length(node.labels[[1]]))
node.labels[[2]] <- rep(19, length(node.labels[[2]]))
node.labels[[3]][grep("astyanax_cave", node.labels[[3]])] = 17
node.labels[[3]][grep("astyanax_surface", node.labels[[3]])] = 19
node.labels[[3]] <- as.numeric(node.labels[[3]])
node.labels[[4]][grep("astyanax_cave", node.labels[[4]])] = 19
node.labels[[4]][grep("astyanax_surface", node.labels[[4]])] = 17
node.labels[[4]] <- as.numeric(node.labels[[4]])
node.labels[[5]][grep("choy_surface", node.labels[[5]])] = 19
node.labels[[5]][grep("molino_cave", node.labels[[5]])] = 17
node.labels[[5]][grep("pachon_cave", node.labels[[5]])] = 17
node.labels[[5]][grep("tinaja_cave", node.labels[[5]])] = 17
node.labels[[5]] <- as.numeric(node.labels[[5]])
node.labels[[6]][grep("choy_surface", node.labels[[6]])] = 19
node.labels[[6]][grep("molino_cave", node.labels[[6]])] = 17
node.labels[[6]][grep("pachon_cave", node.labels[[6]])] = 17
node.labels[[6]][grep("tinaja_cave", node.labels[[6]])] = 17
node.labels[[6]] <- as.numeric(node.labels[[6]])

# Plot dendrograms as list

dendrograms2 <- dendrograms

for (i in 1:length(dendrograms2)) {
	labels(dendrograms2[[i]]) <- sub('.*e ', "", labels(dendrograms2[[i]]))
}

dend.plots <- list()
for (i in 1:length(dendrograms2)) {
	dend.plots[[i]] <- ggplot(dendrograms2[[i]] %>% dendextend::set("labels_col", value = my_colour_palette[[i]]) %>% dendextend::set("labels_cex", value = 0.5) %>% dendextend::set("leaves_pch", node.labels[[i]]) %>% dendextend::set("leaves_col", label.colours[[i]]) %>% dendextend::set("leaves_cex", 1.5), horiz = T, offset_labels = -50) + theme(plot.background = element_rect(fill = "transparent", color = NA))
}


ggsave("Figures/AstMex_subtype_species_dendrogram.png", units = "in", dpi = 200, height = 10, width = 25, limitsize = FALSE)




