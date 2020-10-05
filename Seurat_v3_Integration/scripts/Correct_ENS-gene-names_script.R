library(ape)

grcz10 <- read.gff("~/Downloads/genome_assemblies_genome_gff/ncbi-genomes-2020-08-13/GCF_000002035.5_GRCz10_genomic.gff")
grcz10 <- grcz10[grcz10$type == "gene",]

head(grcz10)
# seqid              source type start   end score strand phase                                                                                                                                                                                                         attributes
# 2   NC_007112.6          BestRefSeq gene  6642 11878    NA      -  <NA>                             ID=gene0;Dbxref=GeneID:192301,ZFIN:ZDB-GENE-020419-25;Name=rpl24;description=ribosomal protein L24;gbkey=Gene;gene=rpl24;gene_biotype=protein_coding;gene_synonym=chunp6908,wu:fa93e03
# 16  NC_007112.6 BestRefSeq%2CGnomon gene 11926 16373    NA      +  <NA>                                ID=gene1;Dbxref=GeneID:386640,ZFIN:ZDB-GENE-031030-11;Name=cep97;description=centrosomal protein 97;gbkey=Gene;gene=cep97;gene_biotype=protein_coding;gene_synonym=lrriq2,zgc:63856
# 75  NC_007112.6              Gnomon gene 18663 23838    NA      +  <NA>                                                                                             ID=gene2;Dbxref=GeneID:100001781,ZFIN:ZDB-GENE-071024-1;Name=nfkbiz;gbkey=Gene;gene=nfkbiz;gene_biotype=protein_coding
# 103 NC_007112.6          BestRefSeq gene 27688 34330    NA      +  <NA>                                 ID=gene3;Dbxref=GeneID:550463,ZFIN:ZDB-GENE-050417-287;Name=eed;description=embryonic ectoderm development;gbkey=Gene;gene=eed;gene_biotype=protein_coding;gene_synonym=zgc:112509
# 129 NC_007112.6          BestRefSeq gene 36734 39191    NA      +  <NA>                                                                ID=gene4;Dbxref=GeneID:550239,ZFIN:ZDB-GENE-050417-34;Name=zgc:110091;description=zgc:110091;gbkey=Gene;gene=zgc:110091;gene_biotype=protein_coding
# 141 NC_007112.6          BestRefSeq gene 39325 44525    NA      -  <NA> ID=gene5;Dbxref=GeneID:323208,ZFIN:ZDB-GENE-030131-1928;Name=tmem39a;description=transmembrane protein 39A;gbkey=Gene;gene=tmem39a;gene_biotype=protein_coding;gene_synonym=fb92e03,wu:fb92e03,zgc:65781,zgc:77444


# Can strsplit the attributes by ";", then by "=", and take the first entry as names?

test <- strsplit(grcz10$attributes, ";")
test2 <- lapply(test, function(x) strsplit(x, "="))
test3 <- list()

for(i in 1:length(test2)) {
  test3[[i]] <- lapply(test2[[i]], function(x) x[[2]])
  names(test3[[i]]) <- lapply(test2[[i]], function(x) x[[1]])
}

test4 <- lapply(test3, function(x) unlist(x))
test5 <- Reduce(rbind, test4)
test6 <- as.data.frame(test5)
test7 <- separate(test6, "Dbxref", sep = ",", into = c("GeneID", "ZFIN")) ## Now need to strsplit by : for GeneID and ZFIN

test7$GeneID <- substr(test7$GeneID, start = 8, stop = max(apply(test7, 1, function(x) nchar(x[2])))+1)
test7$ZFIN <- substr(test7$ZFIN, start = 8, stop = 100)

GRCz10 <- test7

# Can use the above to match NCBI GeneIDs with ensembl ortholog list
mart.ens <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_GRCz10-ENS-NCBI.txt")

mart.zeb <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_GRCz10_orthologs.txt")

mart.zeb.paralogs <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_GRCz10_paralogs.txt")

# nested match and save the output

mart.zeb$GeneID <- GRCz10$Name[match(mart.ens$NCBI.gene.ID[match(mart.zeb$Gene.stable.ID, mart.ens$Gene.stable.ID)], GRCz10$GeneID)]
write.csv(mart.zeb, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_GRCz10_orthologs_corrected.txt")


# Replace Gene.name and Zebrafish.paralog.associated.gene.name with the NCBI ID, if the NCBI ID is in row.names of the Seurat object, but double check that the Gene.name is not in the object
names <- row.names(GetAssayData(hypo.zeb))

# replace paralog names, but not if it's already blank (no paralog)
mart.zeb.paralogs$GeneID <- mart.zeb$GeneID[match(mart.zeb.paralogs$Gene.stable.ID, mart.zeb$Gene.stable.ID)]
mart.zeb.paralogs.2 <- mart.zeb.paralogs

mart.zeb.paralogs.2$Gene.name <- ifelse(mart.zeb.paralogs$Gene.name == "", "", 
                                      ifelse(mart.zeb.paralogs$Gene.name %in% unique(mart.zeb.paralogs$Gene.name)[!(unique(mart.zeb.paralogs$Gene.name) %in% names)], 
                                             ifelse(mart.zeb.paralogs$GeneID == "", mart.zeb.paralogs$Gene.name, mart.zeb.paralogs$GeneID),
                                             mart.zeb.paralogs$Gene.name)
                                      
                                      )


mart.zeb.paralogs.2$Zebrafish.paralogue.associated.gene.name <- ifelse(mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name == "", mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name, 
                                                                     
                                                                     ifelse(mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name %in% unique(mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name)[!(unique(mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name) %in% names)], 
                                                                            ifelse(mart.zeb.paralogs$GeneID == "", mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name, mart.zeb.paralogs$GeneID[match(mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name, mart.zeb.paralogs$Gene.name)]), 
                                                                            mart.zeb.paralogs$Zebrafish.paralogue.associated.gene.name))



write.csv(mart.zeb.paralogs.2, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/mart_export_GRCz10_paralogs_corrected.txt")


