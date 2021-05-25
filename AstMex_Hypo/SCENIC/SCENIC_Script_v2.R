library(reshape2)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(scales)
library(tictoc)
library(ggrepel)

## Re-make for astyanax - no need to load the two objects, but will need to split one of them

## Load arguments

args <- commandArgs(trailingOnly = TRUE)

# For testing script
# setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/")
# args <- c("cave", "synaptic", 500, 1)

if (length(args) != 4) {
	stop("Need four arguments; morph, go list, nTrees, and re/start position")
}
if (args[[1]] != "surface" & args[[1]] != "cave" & args[[1]] != "pachon"  & args[[1]] != "tinaja"  & args[[1]] != "molino") {
	stop("First argument must be 'surface', 'cave', 'pachon', 'tinaja', or 'molino'")
}

# GO_list <- read.csv(paste("/Volumes/BZ/Home/gizevo30/R_Projects/SCENIC_Hypo/", args[[2]], sep = ""), head = F)

if (args[[2]] == "nps") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_neuropeptide.csv", head = F)
}
if (args[[2]] == "synaptic") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_synaptic_vesicle.csv", head = F)
}
if (args[[2]] == "nts") {
	GO_list1 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_glutamine.csv", head = F)
	GO_list2 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_GABA.csv", head = F)
	GO_list3 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_neurotransmitter.csv", head = F)
	GO_list4 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_serotonin.csv", head = F)
	GO_list5 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_dopamine.csv", head = F)
	GO_list6 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_epinephrine.csv", head = F)
	GO_list <- rbind(GO_list1, GO_list2, GO_list3, GO_list4, GO_list5, GO_list6)
	rm(GO_list1, GO_list2, GO_list3, GO_list4, GO_list5, GO_list6)
}
if (args[[2]] == "ion") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_ion_channel.csv", head = F)
}
if (args[[2]] == "behavior") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_behavior.csv", head = F)
}

## Prepare gene lists for Targets and TFs
GO_list <- sort(as.character(unique(GO_list$V2)))
TF_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/GO_TF_list.csv", head = F)
TF_list <- as.character(unique(TF_list$V2))
GO_list <- GO_list[!(GO_list %in% TF_list)]

## load functions
source("SCENIC_functions.R")

## If restarting mid GENIE3, this code chunk does that, then ends the script without overwriting files in /int
args[[4]] <- as.numeric(args[[4]])

if (args[[4]] != 1) {
	setwd(paste("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC", args[[1]], args[[2]], sep = "/"))
	
	if (args[[1]] == "pachon") {
		hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_pachon_5k_vSCENIC.rds")
		exprMat <- as.matrix(GetAssayData(hypo))
		rm(hypo)
	}

	if (args[[1]] == "tinaja") {
	  hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_tinaja_5k_vSCENIC.rds")
		exprMat <- as.matrix(GetAssayData(hypo))
		rm(hypo)
	}
	
  if (args[[1]] == "molino") {
    hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_molino_5k_vSCENIC.rds")
    exprMat <- as.matrix(GetAssayData(hypo))
    rm(hypo)
  }
  
  if (args[[1]] == "surface") {
    hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_surface_20k_vSCENIC.rds")
    exprMat <- as.matrix(GetAssayData(hypo))
    rm(hypo)
  }
  
  if (args[[1]] == "cave") {
    hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_cave_20k_vSCENIC.rds")
    exprMat <- as.matrix(GetAssayData(hypo))
    rm(hypo)
  }
  
	# load previous work
	genesKept <- readRDS("int/1.1_genesKept.Rds")
	# corrMat <- readRDS("int/1.2_corrMat.Rds")
	scenicOptions <- readRDS(file="int/scenicOptions.Rds")
	
	exprMat <- exprMat[genesKept, ]
	
	# recreate weightMatrix file
	parts <- list.files("int/")[grep("part_", list.files("int/"))]
	parts <- paste("int/", parts, sep = "")
	weightMatrices <- lapply(parts, function(x) readRDS(file = x))
	
	# Run GENIE3 on TF list (my own) and the GO list as targets
	GO_list <- GO_list[GO_list %in% rownames(exprMat)]
	
	print(paste("Running Genie3 for", length(GO_list), "across", round(length(GO_list)/8), "parts", sep = " "))
	
	# Run Genie3 and time it's execution (for planning purposes)
	tic(paste("Genie3 with ", nrow(exprMat), " genes, ", ncol(exprMat), " cells, and ", length(TF_list), " TFs", sep = ""))
	runGenie3v2(exprMat, scenicOptions, GO_list = GO_list, nParts = round(length(GO_list)/8), nTrees = as.numeric(args[[3]]), start = args[[4]], weightMatrices = weightMatrices)
	toc()
	
	## Build and score the GRN
	# Run the remaining steps using the wrapper functions:
	scenicOptions <- readRDS("int/scenicOptions.Rds")
	scenicOptions@settings$verbose <- TRUE
	scenicOptions@settings$nCores <- 10
	scenicOptions@settings$seed <- 123
	scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run
	
	tic(paste("SCENIC step 1 with ", nrow(exprMat), " genes, ", ncol(exprMat), " cells, and ", length(TF_list), " TFs", sep = ""))
	runSCENIC_1_coexNetwork2modules(scenicOptions)
	toc()
	
	stop(paste("Finished aftering restarting from GENIE3 step ", args[[4]], sep = ""))
}

## Otherwise make the directories, and start from 1

# Set and create directories
if (!dir.exists(args[[1]])) {
  dir.create(args[[1]])
}

setwd(args[[1]])

if (!dir.exists(args[[2]])) {
	dir.create(args[[2]])
	dir.create(paste(args[[2]], "/int",  sep = ""))
}
setwd(args[[2]])


# Load scenicOptions and other files from previous run, and save to newly directory
# Need to run geneFilteringv2 and runCorrelation and save the output to the nps int folder, where it can be found for subsequent runs 
# nCountsPerGene should be low enough to capture "hcrt", the np we want to keep

# # for surface
# genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=20, minSamples=8)
# saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/surface/nps/int/1.1_genesKept.Rds")
# 
# # for cave
# genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=30, minSamples=11)
# saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cave/nps/int/1.1_genesKept.Rds")
# 
# # For both
# exprMat <- exprMat[genesKept, ]
# 
# tic("correlation")
# runCorrelation(exprMat, scenicOptions)
# toc()


if (args[[1]] == "pachon") {
  hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_pachon_5k_vSCENIC.rds")
  exprMat <- as.matrix(GetAssayData(hypo))
  cellInfo <- data.frame(hypo@meta.data)
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  rm(hypo)
  
  if (args[[2]] == "nps") {
    org="mgi" # or hgnc, or dmel
    dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cisTarget_databases" # RcisTarget databases location
    myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
    
    scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
    saveRDS(scenicOptions, file="int/scenicOptions.Rds")
    
    counts <- rowSums(exprMat, na.rm = T)["hcrt"]
    cells <- rowSums(exprMat > 0, na.rm = T)["hcrt"]
    
    genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=floor(counts-1), minSamples=floor(cells-1))
    saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/surface/nps/int/1.1_genesKept.Rds")
    
    exprMat <- exprMat[genesKept, ]
    
    tic("correlation")
    runCorrelation(exprMat, scenicOptions)
    toc()
  }
  
  genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/pachon/nps/int/1.1_genesKept.Rds")
  saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
  corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/pachon/nps/int/1.2_corrMat.Rds")
  saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}

if (args[[1]] == "tinaja") {
  hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_tinaja_5k_vSCENIC.rds")
  exprMat <- as.matrix(GetAssayData(hypo))
  cellInfo <- data.frame(hypo@meta.data)
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  rm(hypo)
  
  if (args[[2]] == "nps") {
    org="mgi" # or hgnc, or dmel
    dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cisTarget_databases" # RcisTarget databases location
    myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
    
    scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
    saveRDS(scenicOptions, file="int/scenicOptions.Rds")
    
    counts <- rowSums(exprMat, na.rm = T)["hcrt"]
    cells <- rowSums(exprMat > 0, na.rm = T)["hcrt"]
    
    genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=floor(counts-1), minSamples=floor(cells-1))
    saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/surface/nps/int/1.1_genesKept.Rds")
    
    exprMat <- exprMat[genesKept, ]
    
    tic("correlation")
    runCorrelation(exprMat, scenicOptions)
    toc()
  }
  
  genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/tinaja/nps/int/1.1_genesKept.Rds")
  saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
  corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/tinaja/nps/int/1.2_corrMat.Rds")
  saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}

if (args[[1]] == "molino") {
  hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_molino_5k_vSCENIC.rds")
  exprMat <- as.matrix(GetAssayData(hypo))
  cellInfo <- data.frame(hypo@meta.data)
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  rm(hypo)
  
  if (args[[2]] == "nps") {
    org="mgi" # or hgnc, or dmel
    dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cisTarget_databases" # RcisTarget databases location
    myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
    
    scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
    saveRDS(scenicOptions, file="int/scenicOptions.Rds")
    
    counts <- rowSums(exprMat, na.rm = T)["hcrt"]
    cells <- rowSums(exprMat > 0, na.rm = T)["hcrt"]
    
    genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=floor(counts-1), minSamples=floor(cells-1))
    saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cave/nps/int/1.1_genesKept.Rds")
    
    exprMat <- exprMat[genesKept, ]
    
    tic("correlation")
    runCorrelation(exprMat, scenicOptions)
    toc()
  }
  
  genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/molino/nps/int/1.1_genesKept.Rds")
  saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
  corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/molino/nps/int/1.2_corrMat.Rds")
  saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}

if (args[[1]] == "surface") {
  hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_surface_20k_vSCENIC.rds")
  exprMat <- as.matrix(GetAssayData(hypo))
  cellInfo <- data.frame(hypo@meta.data)
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  rm(hypo)
  
  if (args[[2]] == "nps") {
    org="mgi" # or hgnc, or dmel
    dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cisTarget_databases" # RcisTarget databases location
    myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
    
    scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
    saveRDS(scenicOptions, file="int/scenicOptions.Rds")
    
    counts <- rowSums(exprMat, na.rm = T)["hcrt"]
    cells <- rowSums(exprMat > 0, na.rm = T)["hcrt"]
    
    genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=floor(counts-1), minSamples=floor(cells-1))
    saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cave/nps/int/1.1_genesKept.Rds")
    
    exprMat <- exprMat[genesKept, ]
    
    tic("correlation")
    runCorrelation(exprMat, scenicOptions)
    toc()
  }
  
  genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/surface/nps/int/1.1_genesKept.Rds")
  saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
  corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/surface/nps/int/1.2_corrMat.Rds")
  saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}

if (args[[1]] == "cave") {
  hypo <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_cave_20k_vSCENIC.rds")
  exprMat <- as.matrix(GetAssayData(hypo))
  cellInfo <- data.frame(hypo@meta.data)
  saveRDS(cellInfo, file="int/cellInfo.Rds")
  rm(hypo)
  
  if (args[[2]] == "nps") {
    org="mgi" # or hgnc, or dmel
    dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cisTarget_databases" # RcisTarget databases location
    myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
    data(defaultDbNames)
    dbs <- defaultDbNames[[org]]
    scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)
    
    scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
    saveRDS(scenicOptions, file="int/scenicOptions.Rds")
    
    counts <- rowSums(exprMat, na.rm = T)["hcrt"]
    cells <- rowSums(exprMat > 0, na.rm = T)["hcrt"]
    
    genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=floor(counts-1), minSamples=floor(cells-1))
    saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cave/nps/int/1.1_genesKept.Rds")
    
    exprMat <- exprMat[genesKept, ]
    
    tic("correlation")
    runCorrelation(exprMat, scenicOptions)
    toc()
  }
  
  genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cave/nps/int/1.1_genesKept.Rds")
  saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
  corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cave/nps/int/1.2_corrMat.Rds")
  saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}


## Initialize SCENIC settings (using mouse as stand in)

org="mgi" # or hgnc, or dmel
dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/SCENIC/cisTarget_databases" # RcisTarget databases location
myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

## Run the first step of SCENIC, to recover gene correlation networks, without running the downstream filtering and quantifications

## Co-expression network
# Gene filter/selection from loaded files
exprMat <- exprMat[genesKept, ]

# Run GENIE3 on TF list (my own) and the GO list as targets

GO_list <- GO_list[GO_list %in% rownames(exprMat)]

# Run Genie3 and time it's execution (for planning purposes)

warning(paste("Executing Genie3 over", round(length(GO_list)/8), "parts", sep = " "))

tic(paste("Genie3 with ", nrow(exprMat), " genes, ", ncol(exprMat), " cells, and ", length(TF_list), " TFs", sep = ""))
runGenie3v2(exprMat, scenicOptions, GO_list = GO_list, nParts = round(length(GO_list)/8), nTrees = as.numeric(args[[3]]), start = 1)
toc()


## Build and score the GRN
# Run the remaining steps using the wrapper functions:
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123


scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run

tic(paste("SCENIC step 1 with ", nrow(exprMat), " genes, ", ncol(exprMat), " cells, and ", length(TF_list), " TFs", sep = ""))
runSCENIC_1_coexNetwork2modules(scenicOptions)
toc()
