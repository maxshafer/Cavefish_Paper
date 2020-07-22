library("biomaRt")
ensembl <- useMart("ensembl")
ensembl <- useDataset("drerio_gene_ensembl", mart = ensembl) # or tilapia etc "oniloticus_gene_ensembl"

## Retrive zebrafish genes (ensembl ids and gene names) and associated tilapia and stickleback ensembl ids (to match with the data from paper)
mart <- getBM(attributes = c("external_gene_name", "ensembl_gene_id", "mmusculus_homolog_associated_gene_name"), mart = ensembl)


# For loop to make gene lists for Gene Analytics (mouse orthologs of zebrafish genes)

for (i in 0:36) {
	markers <- read.csv(paste("cluster_", i, ".csv", sep = ""))
	markers <- subset(mart, external_gene_name %in% markers$gene)
	write.table(markers$mmusculus_homolog_associated_gene_name, file = paste("CellAnalytics/cluster_", i, ".txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
}


# Subset biomart results by cluster markers
markers <- read.csv("cluster_0.csv")
mart2 <- subset(mart, external_gene_name %in% markers$gene)
write.table(mart2$mmusculus_homolog_associated_gene_name, file = "cluster_0_analytics.txt", row.names = FALSE, quote = FALSE)


for (i in 0:36) {
  markers <- read.csv(paste("cluster_", i, ".csv", sep = ""))
  markers <- as.character(markers$gene)
  write.table(markers, file = paste("CellAnalytics/cluster_", i, ".txt", sep = ""), row.names = FALSE, col.names = FALSE, quote = FALSE)
}
