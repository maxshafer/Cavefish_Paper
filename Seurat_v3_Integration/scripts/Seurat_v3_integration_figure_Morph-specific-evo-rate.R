library(alluvial)
library(Seurat)
library(reshape)
library(scales)
library(networkD3)

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/")


df <- data.frame(branch = c("zebrafish", "LCA A. mex", "LCA Molino + Surface", "Surface", "Molino", "LCA Pachon + Tinaja", "Pachon", "Tinaja"), 
                 cell_type_number = c(12, 25, 0, 0, 33, 4, 8, 13),
                 time_upper = c(150000000, 149500000, 340000, 160000, 160000, 270000, 230000, 230000))

df$rate_per100ky <- df$cell_type_number/df$time_upper*100000
df$log_rate_per100ky <- log(df$rate_per100ky)

number <- ggplot(df, aes(x = branch, y = cell_type_number, colour = cell_type_number)) + geom_point(size = 5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_colour_viridis()
per100k <- ggplot(df, aes(x = branch, y = log_rate_per100ky, colour = log_rate_per100ky)) + geom_point(size = 5) + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) + scale_colour_viridis()


dev.new()
number + per100k

## Use these plots to colour the branchs of the cladogram in Affinity