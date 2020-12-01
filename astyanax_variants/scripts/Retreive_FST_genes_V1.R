library(Seurat)
library(data.table)
library(Gviz)
library(biomaRt)
library(qqman)
library(DescTools)
library(stringr)
library(fdrtool)
library(dplyr)
library(data.table)
library(purrr)
library(pbapply)
library(venneuler)
# library(refGenome) #


########################################################
# Generate lists of dataframes for each
########################################################

setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/astyanax_variants/FST")
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/astyanax_variants/scripts/Astyanax_variants_functions.R")

library(Rgb)
gtf <- read.gtf("../Astyanax_mexicanus.AstMex102.91.gtf")
gtf2 <- gtf[gtf$feature == "gene",]


# Retrieve dataset names for weir.fst
datasets.weir <- list.files()[grep("windowed.weir",list.files())]
# datasets.weir <- list.files()[grep("perbase",list.files())]

# Generate weir.fst tables from datasets list
weir.fst <- lapply(datasets.weir, readFST)

# Determine and apply fdr cutoff (0.5)
weir.fst.top <- lapply(weir.fst, function(x) fdrcutoff(x, cutoff = 0.05))
weir.fst.top.01 <- lapply(weir.fst, function(x) fdrcutoff(x, cutoff = 0.001))

# Take only those with outlier high values (not low)
weir.fst.top[1:3] <- lapply(weir.fst.top[1:3], function(x) x[x$MEAN_FST > 0.1, ])
weir.fst.top[4:6] <- lapply(weir.fst.top[4:6], function(x) x[x$MEAN_FST > 0.2, ])
weir.fst.top.01 [1:3] <- lapply(weir.fst.top.01 [1:3], function(x) x[x$MEAN_FST > 0.1, ])
weir.fst.top.01 [4:6] <- lapply(weir.fst.top.01 [4:6], function(x) x[x$MEAN_FST > 0.2, ])

# ## Alternatively, take only the top % of FST values
# weir.fst.top.2 <- lapply(weir.fst, function(x) tibble(x %>% top_n(round(.01*nrow(x)), FisherZMean)))

weir.genes <- pblapply(weir.fst.top, function(x) bind_rows(apply(x, 1, getGenesEns)))
weir.genes.01 <- pblapply(weir.fst.top.01, function(x) bind_rows(apply(x, 1, getGenesEns)))

weir.genes.all <- pblapply(weir.fst, function(x) bind_rows(apply(x, 1, getGenesEns)))

# getGenesEns introduces NAs for gene_ids with no gene_name, convert to gene_ids
for (i in 1:length(weir.genes)) {
  weir.genes[[i]]$gene_name[is.na(weir.genes[[i]]$gene_name)] <- weir.genes[[i]]$gene_id[is.na(weir.genes[[i]]$gene_name)]
}

for (i in 1:length(weir.genes.all)) {
  weir.genes.all[[i]]$gene_name[is.na(weir.genes.all[[i]]$gene_name)] <- weir.genes.all[[i]]$gene_id[is.na(weir.genes.all[[i]]$gene_name)]
}

names(weir.genes) <- datasets.weir
names(weir.genes.all) <- datasets.weir
saveRDS(weir.genes, file = "../weir.genes.rds")
saveRDS(weir.genes.all, file = "../weir.genes.all.rds")

weir.genes <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/astyanax_variants/weir.genes.rds")



##### MAKE PLOTS #######

## Find union of genes associated with divergent FST windows between Pach/Tin and surface (INDELs and SNPS)
snp.genes <- Reduce(union, lapply(weir.genes, function(x) x$gene_name))

snp.genes.molino <- union(weir.genes[[1]]$gene_name, weir.genes[[4]]$gene_name)
snp.genes.pachon <- union(weir.genes[[2]]$gene_name, weir.genes[[5]]$gene_name)
snp.genes.tinaja <- union(weir.genes[[3]]$gene_name, weir.genes[[6]]$gene_name)

## Load marker, sub markers, and differentially expressed genes from file
gene.lists <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/marker_gene_lists.rds")

marker.genes <- Reduce(union, lapply(gene.lists$conserved.markers.sub, function(x) row.names(x)))
de.genes1 <- lapply(gene.lists[[8]], function(x) row.names(x[x$p_val_adj < 0.05,]))
de.genes <- Reduce(union, lapply(gene.lists[[8]], function(x) row.names(x[x$p_val_adj < 0.05,])))
de.genes2 <- lapply(gene.lists[[8]], function(x) x[row.names(x) %in% snp.genes,])

overlap.genes <- union(snp.genes[snp.genes %in% marker.genes], snp.genes[snp.genes %in% de.genes])

list <- list(markers = marker.genes, de = de.genes, snp = snp.genes, snp.pach = snp.genes.pachon, snp.tinaja = snp.genes.tinaja, snp.molino = snp.genes.molino)

# Make venn diagrams

venn.data <- c(Marker_genes = length(list$markers[!(list$markers %in% c(list$snp, list$de))]), 
               DE_genes = length(list$de[!(list$de %in% c(list$markers, list$snp))]), 
               Fst_genes = length(list$snp[!(list$snp %in% c(list$markers, list$de))]), 
               'Marker_genes&DE_genes' = length(intersect(list$markers, list$de)[!(intersect(list$markers, list$de) %in% list$snp)]),
               'Marker_genes&Fst_genes' = length(intersect(list$markers, list$snp)[!(intersect(list$markers, list$snp) %in% list$de)]),
               'DE_genes&Fst_genes' = length(intersect(list$de, list$snp)[!(intersect(list$de, list$snp) %in% list$markers)]),
               'Marker_genes&DE_genes&Fst_genes' = length(Reduce(intersect, list(list$markers, list$de, list$snp))))

vd1 <- venneuler(venn.data)

plot(vd1)

venn.data.morphs <- c(Pachon = length(list$snp.pach[!(list$snp.pach %in% c(list$snp.tinaja, list$snp.molino))]),
                      Tinaja = length(list$snp.tinaja[!(list$snp.tinaja %in% c(list$snp.pach, list$snp.molino))]),
                      Molino = length(list$snp.molino[!(list$snp.molino %in% c(list$snp.pach, list$snp.tinaja))]),
                      'Pachon&Tinaja' = length(intersect(list$snp.pach, list$snp.tinaja)[!(intersect(list$snp.pach, list$snp.tinaja) %in% list$snp.molino)]),
                      'Pachon&Molino' = length(intersect(list$snp.pach, list$snp.molino)[!(intersect(list$snp.pach, list$snp.molino) %in% list$snp.tinaja)]),
                      'Tinaja&Molino' = length(intersect(list$snp.molino, list$snp.tinaja)[!(intersect(list$snp.molino, list$snp.tinaja) %in% list$snp.pach)]),
                      'Pachon&Tinaja&Molino' = length(Reduce(intersect, list(list$snp.pach, list$snp.molino, list$snp.tinaja))))

vd2 <- venneuler(venn.data.morphs)
plot(vd2)


david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")
setAnnotationCategories(david, c("UP_KEYWORDS", "GOTERM_BP_DIRECT", "GOTERM_CC_DIRECT", "GOTERM_MF_DIRECT", "KEGG_PATHWAY"))

david.genes <- select(org.Dr.eg.db, overlap.genes, "ENTREZID", "SYMBOL")[,2]
result <- addList(david, david.genes, idType = "ENTREZ_GENE_ID", listName = "Glut_1_4", listType = "Gene")
annoCharts <- getFunctionalAnnotationChart(david)

anno.data <- annoCharts[,c(1,2,3,5,7,8,9,10,11,12,13)]

# Make DAVID plot

david.plot <- ggplot(anno.data[anno.data$PValue < 0.01,], aes(x=Term, y=Count, fill=log(Bonferroni))) + geom_bar(stat="identity") + coord_flip() + theme_classic() + theme(axis.text.y = element_blank())

david.plot <- david.plot + geom_text(aes(x=Term, y=Count+0.5, label = Term), hjust = 0) + theme(legend.position = c(0.75,0.75), legend.background = element_blank())
david.plot



library(formattable)

# Get Go lists and use only those that are marker genes

go_lists <- list.files("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/")[grep("GO", list.files("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/"))]
go_lists <- lapply(go_lists, function(x) read.csv(paste("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/", x, sep = ""), head = F))
go_lists <- lapply(go_lists, function(x) unique(x$V2))

names(go_lists) <- list.files("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/")[grep("GO", list.files("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/GO_lists/"))]
names(go_lists)

circadian.genes <- c("slc1a1b", "galn", go_lists$GO_circadian.csv)
np.genes <- c("trh", go_lists$GO_neuropeptide.csv)
genes <- c(circadian.genes, np.genes)

df <- list()
for (i in 1:length(genes)) {
  index <- lapply(weir.genes, function(x) grep(genes[[i]], x$gene_name))
  df[[i]] <- Reduce(rbind, (lapply(seq_along(index), function(x) {
    if(length(index[[x]]) > 0) {
      dft <- weir.genes[[x]][index[[x]],]
      dft$comparison <- names(weir.genes)[[x]]
      return(dft)
    }
  })))
}

df <- Reduce(rbind, df)
df$Type[grep("snps", df$comparison)] <- "SNPs"
df$Type[grep("indels", df$comparison)] <- "INDELs"
df$Comparison[grep("tinaja", df$comparison)] <- "Tinaja"
df$Comparison[grep("pachon", df$comparison)] <- "Pachon"
df$Comparison[grep("molino", df$comparison)] <- "Molino"
df <- df[df$gene_name %in% genes, c("gene_name", "Comparison", "Type", "CHROM", "BIN_START", "BIN_END", "MEAN_FST", "N_VARIANTS")]

# Take rows for Maximum MEAN_FST based on columns c(1,2,3)

df2 <- df %>% group_by(gene_name, Comparison, Type, CHROM) %>% dplyr::summarize(BIN_START = min(BIN_START), BIN_END = max(BIN_END), MEAN_FST = mean(MEAN_FST), N_VARIANTS = sum(N_VARIANTS))

df3 <- dcast(setDT(df2), gene_name+CHROM+BIN_START+BIN_END~Comparison+Type, value.var=c("MEAN_FST", "N_VARIANTS"))
df3[is.na(df3)] <- 0
df3$N_Variants <- apply(df3, 1, function(x) sum(as.numeric(x[c(11:16)])))
df3$GO_Term <- ifelse(df3$gene_name %in% circadian.genes & df3$gene_name %in% np.genes, "Circadian/NP", ifelse(df3$gene_name %in% circadian.genes, "Circadian", "Neuropeptide"))

df3 <- df3[,c(1,18,2:4,17,5:10)]

df3[,6:11] <- round(df3[,6:11], digits = 3)
df3[df3 == 0] <- ""
colnames(df3) <- c("Gene", "GO_Term", "Chromosome", "Bin start", "Bin end", "Variants", "M_I", "M_S", "P_I", "P_S", "T_I", "T_S")

df4 <-formattable(df3, list(
  M_I = color_tile("white", "orange"),
  M_S = color_tile("white", "orange"),
  P_I = color_tile("white", "orange"),
  P_S = color_tile("white", "orange"),
  T_I = color_tile("white", "orange"),
  T_S = color_tile("white", "orange"),
  GO_Term = color_tile("#FDE725FF", "#22A884FF"),
  area(col = c(Variants)) ~ normalize_bar("pink", 0.2)
))

# Need to save as html (from R), then use website to convert to PDF (https://pdfcrowd.com/#convert_by_upload+with_options)

df4









# Should make a list of the DE genes between cave/surface for each subcluster

# for each subcluster, the surface and cave specific genes
test <- lapply(names(gene.lists$conserved.markers.sub), function(x) union(setdiff(row.names(gene.lists$surface.markers.sub[[x]]), row.names(gene.lists$conserved.markers.sub[[x]])), setdiff(row.names(gene.lists$cave.markers.sub[[x]]), row.names(gene.lists$conserved.markers.sub[[x]]))))
names(test) <- names(gene.lists$conserved.markers.sub)
test2 <- lapply(test, function(x) x[x %in% snp.genes])

unlist(lapply(seq_along(test), function(x) length(test2[[x]])/length(test[[x]])))

combined.markers <- lapply(names(gene.lists$conserved.markers.sub), function(x) unlist(Reduce(union, list(row.names(gene.lists$conserved.markers.sub[[x]]), row.names(gene.lists$cave.markers.sub[[x]]), row.names(gene.lists$surface.markers.sub[[x]])))))
snp.markers <- lapply(seq_along(gene.lists$conserved.markers.sub), function(x) snp.genes[snp.genes %in% combined.markers[[x]]])

de.markers <- lapply(names(gene.lists$conserved.markers.sub), function(x) de.genes1[[x]])
snp.de <- lapply(seq_along(gene.lists$conserved.markers.sub), function(x) snp.genes[snp.genes %in% de.markers[[x]]])



per.snp.markers <- unlist(lapply(snp.markers, function(x) length(x)))/unlist(lapply(combined.markers, function(x) length(x)))

snp.df <- data.frame(snp.genes = unlist(lapply(snp.markers, function(x) length(x))), marker.genes = unlist(lapply(combined.markers, function(x) length(x))), de.genes = unlist(lapply(de.markers, function(x) length(x))), snp.de = unlist(lapply(snp.de, function(x) length(x))))
snp.df$percent.snp <- snp.df$snp.genes/snp.df$marker.genes*100
snp.df$de.percent.snp <- snp.df$snp.de/snp.df$de.genes*100
snp.df$cell_type <- names(gene.lists$conserved.markers.sub)

ggplot(reshape2::melt(snp.df), aes(x = cell_type, y = value, group = variable, color = variable)) + geom_point() + geom_line() + theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6), axis.title = element_blank()) + facet_wrap(~variable,scales = "free_y", nrow = 6)


