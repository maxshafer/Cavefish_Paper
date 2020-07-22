library(RDAVIDWebService)
library(org.Dr.eg.db)

# Create the david object, associated with registered email address
david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")

# Examine object, and any gene lists associated with it
david

# Potential ID types for gene lists

getIdTypes(david)
# [1] "AFFYMETRIX_3PRIME_IVT_ID" "AFFYMETRIX_EXON_ID"      
# [3] "AGILENT_CHIP_ID"          "AGILENT_ID"              
# [5] "AGILENT_OLIGO_ID"         "APHIDBASE_ID"            
# [7] "BEEBASE_ID"               "BEETLEBASE_ID"           
# [9] "BGD_ID"                   "CGNC_ID"                 
# [11] "CRYPTODB_ID"              "DICTYBASE_ID"            
# [13] "ENSEMBL_GENE_ID"          "ENSEMBL_TRANSCRIPT_ID"   
# [15] "ENTREZ_GENE_ID"           "FLYBASE_GENE_ID"         
# [17] "GENBANK_ACCESSION"        "GENOMIC_GI_ACCESSION"    
# [19] "GENPEPT_ACCESSION"        "LOCUS_TAG"               
# [21] "MGI_ID"                   "MIRBASE_ID"              
# [23] "MRNA_GI_ACCESSION"        "NASONIABASE_ID"          
# [25] "PROTEIN_GI_ACCESSION"     "PSEUDOCAP_ID"            
# [27] "REFSEQ_MRNA"              "REFSEQ_PROTEIN"          
# [29] "RGD_ID"                   "SGD_ID"                  
# [31] "TAIR_ID"                  "UNIGENE"                 
# [33] "UNIPROT_ACCESSION"        "UNIPROT_ID"              
# [35] "VECTORBASE_ID"            "WORMBASE_GENE_ID"        
# [37] "XENBASE_ID"               "ZFIN_ID"

# Annotation categories on DAVID

getAllAnnotationCategoryNames(david)
# [1] "ENTREZ_GENE_ID"            "BBID"                     
# [3] "BIOCARTA"                  "BIOGRID_INTERACTION"      
# [5] "CGAP_EST_QUARTILE"         "CGAP_SAGE_QUARTILE"       
# [7] "CHROMOSOME"                "COG_ONTOLOGY"             
# [9] "CYTOBAND"                  "DIP"                      
# [11] "ENSEMBL_GENE_ID"           "EC_NUMBER"                
# [13] "GAD_DISEASE"               "ENTREZ_GENE_SUMMARY"      
# [15] "GENE3D"                    "GAD_DISEASE_CLASS"        
# [17] "GNF_U133A_QUARTILE"        "GENERIF_SUMMARY"          
# [19] "GOTERM_BP_2"               "GOTERM_BP_1"              
# [21] "GOTERM_BP_4"               "GOTERM_BP_3"              
# [23] "GOTERM_BP_ALL"             "GOTERM_BP_5"              
# [25] "GOTERM_BP_FAT"             "GOTERM_BP_DIRECT"         
# [27] "GOTERM_CC_3"               "GOTERM_CC_4"              
# [29] "GOTERM_CC_1"               "GOTERM_CC_2"              
# [31] "GOTERM_CC_DIRECT"          "GOTERM_CC_FAT"            
# [33] "GOTERM_CC_5"               "GOTERM_CC_ALL"            
# [35] "GOTERM_MF_3"               "GOTERM_MF_4"              
# [37] "GOTERM_MF_1"               "GOTERM_MF_2"              
# [39] "GOTERM_MF_DIRECT"          "GOTERM_MF_FAT"            
# [41] "GOTERM_MF_5"               "GOTERM_MF_ALL"            
# [43] "HIV_INTERACTION_PUBMED_ID" "HIV_INTERACTION_CATEGORY" 
# [45] "HIV_INTERACTION"           "KEGG_PATHWAY"             
# [47] "INTERPRO"                  "INTACT"                   
# [49] "OMIM_DISEASE"              "MINT"                     
# [51] "PIR_SEQ_FEATURE"           "PIR_SUMMARY"              
# [53] "PIR_SUPERFAMILY"           "PFAM"                     
# [55] "PUBMED_ID"                 "REACTOME_PATHWAY"         
# [57] "SMART"                     "PRINTS"                   
# [59] "PRODOM"                    "PROSITE"                  
# [61] "UCSC_TFBS"                 "UP_KEYWORDS"              
# [63] "UNIGENE_EST_QUARTILE"      "SP_COMMENT_TYPE"          
# [65] "SP_COMMENT"                "TIGRFAMS"                 
# [67] "SUPFAM"                    "UP_TISSUE"                
# [69] "UP_SEQ_FEATURE" 

setAnnotationCategories(david, c("GOTERM_BP_ALL","GOTERM_MF_ALL", "GOTERM_CC_ALL","UP_KEYWORDS","KEGG_PATHWAY"))
termCluster <- getClusterReport(david, type="Term")
getClusterReportFile(david, type="Term", fileName="termClusterReport1.tab")
termCluster
head(summary(termCluster))
clustNumber <- 2

# Plot the results of clustering, 1 cluster at a time (genes by go terms)
plot2D(termCluster, clustNumber)


############################################################
############################################################
############################################################


# For loop for making annotation files! Works! Though problematic when using the same david object

david <- DAVIDWebService(email="maxshafer@fas.harvard.edu", url = "https://david.ncifcrf.gov/webservice/services/DAVIDWebService")

# setAnnotationCategories(david, c("UP_KEYWORDS","KEGG_PATHWAY", "COG_ONTOLOGY"))
setAnnotationCategories(david, c("GOTERM_BP_ALL","GOTERM_MF_ALL", "GOTERM_CC_ALL"))


for (i in 0:36) {
  markers <- read.csv(paste("cluster_", i, ".csv", sep = ""))
  markers <- as.character(markers$gene)
  markers <- select(org.Dr.eg.db, markers, "ENTREZID", "SYMBOL")[,2]
  result <- addList(david, markers, idType = "ENTREZ_GENE_ID", listName = paste("cluster", i, sep = ""), listType = "Gene")
  annoChart <- getFunctionalAnnotationChart(david)
  getFunctionalAnnotationChartFile(david, paste("Annotations/cluster_", i, "_annotation_table.csv", sep = ""))
  pdf(paste("Figures/GOTERM/Annotations/cluster_", i, "_annotation.pdf", sep = ""))
  print(plot2D(annoChart, color = c("FALSE" = "black", "TRUE" = "green")))
  dev.off()
  termCluster <- getClusterReport(david, type="Term")
  getClusterReportFile(david, type="Term", fileName="termClusterReport1.tab")
  pdf(paste("Figures/GOTERM/termCluster/cluster_", i, "_termCluster1.pdf", sep = ""))
  print(plot2D(termCluster, 1))
  dev.off()
  pdf(paste("Figures/GOTERM/termCluster/cluster_", i, "_termCluster2.pdf", sep = ""))
  print(plot2D(termCluster, 2))
  dev.off()
  pdf(paste("Figures/GOTERM/termCluster/cluster_", i, "_termCluster3.pdf", sep = ""))
  print(plot2D(termCluster, 3))
  dev.off()
}

plot2D(annoChart, color=c("FALSE"="black", "TRUE"="green"))

markers <- as.character((read.csv(paste("cluster_", i, ".csv", sep = "")))$gene)
markers <- as.character((read.csv("cluster_1.csv"))$gene)


# Load in gene list(s)
cluster_19 <- read.csv("cluster_19.csv")

# The 'gene' column has the info needed
cluster_19 <- as.character(cluster_19$gene)

# Retrieve Entrezid using the symbol as query from the annotated zebrafish genome

library("org.Dr.eg.db")

cluster_19 <- select(org.Dr.eg.db, cluster_19, "ENTREZID", "SYMBOL")[,2]

# Submit as list to DAVID

result <- addList(david, cluster_19, idType = "ENTREZ_GENE_ID", listName = "cluster_19", listType = "Gene")

# Table of genes and their functional annotation
annoTable <- getFunctionalAnnotationTable(david)
getFunctionalAnnotationTableFile(david, "FunctionalAnnotationTable")

## Functional annotation chart (regular type output from DAVID online), with category, term, # genes (count),
## X., PValue, Genes, List.Total, Pop.Hits, Pop.Total, Fold.Enrichment,   Bonferroni,    Benjamini

annoChart <- getFunctionalAnnotationChart(david)
head(annoChart)

annoChart[1:10,1:5]
# Category                 Term Count        X.        PValue
# 1   UP_KEYWORDS    Ribosomal protein    70 23.333333 9.156032e-100
# 2  KEGG_PATHWAY    dre03010:Ribosome    79 26.333333  2.376062e-89
# 3   UP_KEYWORDS    Ribonucleoprotein    70 23.333333  9.098273e-83
# 4  KEGG_PATHWAY    dre04142:Lysosome    18  6.000000  1.786243e-06
# 5   UP_KEYWORDS    Elongation factor     6  2.000000  1.041200e-05
# 6   UP_KEYWORDS      Stress response     7  2.333333  1.098224e-05
# 7  KEGG_PATHWAY   dre04145:Phagosome    14  4.666667  9.806047e-04
# 8   UP_KEYWORDS Protein biosynthesis     7  2.333333  2.166258e-03
# 9   UP_KEYWORDS            Cytoplasm    24  8.000000  3.516369e-03
# 10  UP_KEYWORDS         rRNA-binding     3  1.000000  4.558359e-03



# Gene list, with full names, species 
genelistreport <- getGeneListReport(david)






