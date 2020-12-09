library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)


# Load ense89 mart orthologs

mart.zeb.89 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/ens89_mart_export_GRCz10_orthologs.txt")

# Load corrected ortholog list from SCORPiOS
corr.orthologs <- read.csv("~/Documents/Python_projects/scorpios-workflow/SCORPIOS/SCORPiOs_ensembl89_final/Synteny/Sorted_SyntenyOrthoPred_Clupeocephala_Lepisosteus.oculatus_0", sep = "\t", header = F)
colnames(corr.orthologs) <- c("Cave.fish.gene.stable.ID", "Zebrafish.gene.stable.ID", "synteny_score", "outgroup_gene")

# Filter corrected orthologs for only A.mex + D. rerio orthologs
index <- grep("Astyanax.mexicanus", corr.orthologs$Cave.fish.gene.stable.ID)
orthologs <- corr.orthologs[index,]

# Duplicate lines with mulitple entries in either column 1 or 2

orthologs <- orthologs %>% mutate(Zebrafish.gene.stable.ID = strsplit(as.character(Zebrafish.gene.stable.ID), ",")) %>% unnest(Zebrafish.gene.stable.ID)

orthologs <- orthologs %>% mutate(Cave.fish.gene.stable.ID = strsplit(as.character(Cave.fish.gene.stable.ID), ",")) %>% unnest(Cave.fish.gene.stable.ID)

# Remove species tags
orthologs$Cave.fish.gene.stable.ID <- str_sub(orthologs$Cave.fish.gene.stable.ID, end = -20)
orthologs$Zebrafish.gene.stable.ID <- str_sub(orthologs$Zebrafish.gene.stable.ID, end = -13)

# Replace "None" entries with blanks
orthologs$Cave.fish.gene.stable.ID[orthologs$Cave.fish.gene.stable.ID == "None"] <- ""
orthologs$Zebrafish.gene.stable.ID[orthologs$Zebrafish.gene.stable.ID == "None"] <- ""

orthologs$Cave.fish.gene.name <- mart.zeb.89$Cave.fish.gene.name[match(orthologs$Cave.fish.gene.stable.ID, mart.zeb.89$Cave.fish.gene.stable.ID)]
orthologs$Zebrafish.gene.name <- mart.zeb.89$Gene.name[match(orthologs$Zebrafish.gene.stable.ID, mart.zeb.89$Gene.stable.ID)]

# Add ENS ID if no Cave.fish.gene.name
orthologs$Cave.fish.gene.name[is.na(orthologs$Cave.fish.gene.name)] <- ""
orthologs$Zebrafish.gene.name[is.na(orthologs$Zebrafish.gene.name)] <- ""

index <- orthologs$Cave.fish.gene.name == ""
index2 <- c(1:nrow(orthologs))[index]

orthologs$Cave.fish.gene.name[index2] <- orthologs$Cave.fish.gene.stable.ID[index2]


# Compare the orthologs and mart.zeb.89 dfs

mart.89 <- mart.zeb.89[,c(3,1,4,2)]
colnames(mart.89) <- c("Cave.fish.gene.stable.ID", "Zebrafish.gene.stable.ID", "Cave.fish.gene.name", "Zebrafish.gene.name")

# mart.89 <- mart.89[mart.89$Cave.fish.gene.stable.ID %in% orthologs$Cave.fish.gene.stable.ID,]


library(plyr)
orthologs.2 <- plyr::match_df(orthologs, mart.89, on = c("Cave.fish.gene.stable.ID"))

orthologs.2$ens.89 <- mart.89$Zebrafish.gene.stable.ID[match(orthologs.2$Cave.fish.gene.stable.ID, mart.89$Cave.fish.gene.stable.ID)]

orthologs.2 <- orthologs.2[orthologs.2$Cave.fish.gene.stable.ID != "",]
orthologs.2 <- orthologs.2[orthologs.2$Zebrafish.gene.stable.ID != "",]

orthologs.2$match <- ifelse(orthologs.2$Zebrafish.gene.stable.ID == orthologs.2$ens.89, "yes", "no")

orthologs.2$match2 <- ifelse(orthologs.2$Zebrafish.gene.name == orthologs.2$Cave.fish.gene.name, "yes", "no")

orthologs.3 <- orthologs.2[orthologs.2$match == "no",]
orthologs.3 <- orthologs.3[orthologs.3$match2 == "no",]

orthologs.list <- list()
orthologs.list[[1]] <- orthologs
orthologs.list[[2]] <- orthologs.2
orthologs.list[[3]] <- orthologs.3

saveRDS(orthologs.list, "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/synteny_orthologs.rds")


## test whether zeb or ast paralogs are in the list

zeb <- lapply(gene.lists.para$zebrafish.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.3$Zebrafish.gene.name])
zeb.all <- lapply(gene.lists.para$zebrafish.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.2$Zebrafish.gene.name])

ast <- lapply(gene.lists.para$astyanax.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.3$Cave.fish.gene.name])
ast.all <- lapply(gene.lists.para$astyanax.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.2$Cave.fish.gene.name])


test.df <- data.frame(zeb.para = unlist(lapply(gene.lists.para$zebrafish.markers.sub, function(x) length(x))), 
                      zeb.all = unlist(lapply(zeb.all, function(x) length(x))),
                      zeb.corr = unlist(lapply(zeb, function(x) length(x))), 
                      ast.para = unlist(lapply(gene.lists.para$astyanax.markers.sub, function(x) length(x))),
                      ast.all = unlist(lapply(ast.all, function(x) length(x))),
                      ast.corr = unlist(lapply(ast, function(x) length(x))))

test.df$zeb.percent <- test.df$zeb.corr/test.df$zeb.all
test.df$ast.percent <- test.df$ast.corr/test.df$ast.all


## Test whether zeb or ast markers are in the list

zeb <- lapply(gene.lists.pos$zebrafish.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.3$Zebrafish.gene.name])
zeb.all <- lapply(gene.lists.pos$zebrafish.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.2$Zebrafish.gene.name])

ast <- lapply(gene.lists.pos$astyanax.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.3$Cave.fish.gene.name])
ast.all <- lapply(gene.lists.pos$astyanax.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.2$Cave.fish.gene.name])

con <- lapply(gene.lists.pos$conserved.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.3$Cave.fish.gene.name])
con.all <- lapply(gene.lists.pos$conserved.markers.sub, function(x) row.names(x)[row.names(x) %in% orthologs.2$Cave.fish.gene.name])

test.df.2 <- data.frame(zeb.para = unlist(lapply(gene.lists.pos$zebrafish.markers.sub, function(x) length(x))), 
                      zeb.all = unlist(lapply(zeb.all, function(x) length(x))),
                      zeb.corr = unlist(lapply(zeb, function(x) length(x))), 
                      ast.para = unlist(lapply(gene.lists.pos$astyanax.markers.sub, function(x) length(x))),
                      ast.all = unlist(lapply(ast.all, function(x) length(x))),
                      ast.corr = unlist(lapply(ast, function(x) length(x))),
                      con.all = unlist(lapply(con.all, function(x) length(x))),
                      con.corr = unlist(lapply(con, function(x) length(x))))

test.df.2$zeb.percent <- test.df.2$zeb.corr/test.df.2$zeb.all
test.df.2$ast.percent <- test.df.2$ast.corr/test.df.2$ast.all
test.df.2$con.percent <- test.df.2$con.corr/test.df.2$con.all


data.frame <- data.frame(con.markers = test.df.2$con.percent, zeb.markers = test.df.2$zeb.percent, ast.markers = test.df.2$ast.percent, zeb.para.markers = test.df$zeb.percent, ast.para.markers = test.df$ast.percent)
data.frame$subtype <- names(gene.lists.pos$conserved.markers.sub)

saveRDS(data.frame, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/synteny-corrected-percentage.rds")

