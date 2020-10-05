library(reshape2)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)
library(scales)
library(tictoc)

## Load arguments

args <- commandArgs(trailingOnly = TRUE)

# # For testing script
# setwd("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/")
# args <- c("drerio", "nps", 500, 1)

if (length(args) != 4) {
	stop("Need four arguments; species, go list, nTrees, and re/start position")
}
if (args[[1]] != "drerio" & args[[1]] != "amexicanus") {
	stop("First argument must be 'drerio' or 'amexicanus'")
}

# GO_list <- read.csv(paste("/Volumes/BZ/Home/gizevo30/R_Projects/SCENIC_Hypo/", args[[2]], sep = ""), head = F)

if (args[[2]] == "nps") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_neuropeptide.csv", head = F)
}
if (args[[2]] == "synaptic") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_synaptic_vesicle.csv", head = F)
}
if (args[[2]] == "nts") {
	GO_list1 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_glutamine.csv", head = F)
	GO_list2 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_GABA.csv", head = F)
	GO_list3 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_neurotransmitter.csv", head = F)
	GO_list4 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_serotonin.csv", head = F)
	GO_list5 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_dopamine.csv", head = F)
	GO_list6 <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_epinephrine.csv", head = F)
	GO_list <- rbind(GO_list1, GO_list2, GO_list3, GO_list4, GO_list5, GO_list6)
	rm(GO_list1, GO_list2, GO_list3, GO_list4, GO_list5, GO_list6)
}
if (args[[2]] == "ion") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_ion_channel.csv", head = F)
}
if (args[[2]] == "behavior") {
	GO_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_behavior.csv", head = F)
}

## Prepare gene lists for Targets and TFs
GO_list <- as.character(unique(GO_list$V2))
TF_list <- read.csv("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/GO_DNA_binding.csv", head = F)
TF_list <- as.character(unique(TF_list$V2))

## load functions
source("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/SCENIC_functions.R")

## If restarting mid GENIE3, this code chunk does that, then ends the script without overwriting files in /int
args[[4]] <- as.numeric(args[[4]])

if (args[[4]] != 1) {
	setwd(paste("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC", args[[1]], args[[2]], sep = "/"))
	
	if (args[[1]] == "drerio") {
	  hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
		exprMat <- as.matrix(GetAssayData(hypo.zeb))
		rm(hypo.zeb)
	}

	if (args[[1]] == "amexicanus") {
	  hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")
		exprMat <- as.matrix(GetAssayData(hypo.ast))
		rm(hypo.ast)
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
dir.create(paste("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/", args[[1]], sep = ""))

setwd(paste("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/", args[[1]], sep = ""))

if (!dir.exists(args[[2]])) {
	dir.create(args[[2]])
	dir.create(paste(args[[2]], "/int",  sep = ""))
}
setwd(args[[2]])

# Load scenicOptions and other files from previous run, and save to newly directory
# Need to run geneFilteringv2 and runCorrelation and save the output to the nps int folder, where it can be found for subsequent runs
# nCountsPerGene should be low enough to capture "hcrt", the np we want to keep

# Load scenicOptions and other files from previous run, and save to newly directory

if (args[[1]] == "drerio") {
  hypo.zeb <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/DanRer_Hypo/DanRer_65k.rds")
	exprMat <- as.matrix(GetAssayData(hypo.zeb))
	cellInfo <- data.frame(hypo.zeb@meta.data)
	saveRDS(cellInfo, file="int/cellInfo.Rds")
	rm(hypo.zeb)
	
	if (args[[2]] == "nps") {
	  org="mgi" # or hgnc, or dmel
	  dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/cisTarget_databases" # RcisTarget databases location
	  myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
	  data(defaultDbNames)
	  dbs <- defaultDbNames[[org]]
	  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

	  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
	  saveRDS(scenicOptions, file="int/scenicOptions.Rds")

	  counts <- rowSums(exprMat, na.rm = T)["hcrt"]
	  cells <- rowSums(exprMat > 0, na.rm = T)["hcrt"]

	  genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=floor(counts-1), minSamples=floor(cells-1))
	  saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/drerio/nps/int/1.1_genesKept.Rds")

	  exprMat <- exprMat[genesKept, ]

	  tic("correlation")
	  runCorrelation(exprMat, scenicOptions)
	  toc()
	}
	
	genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/drerio/nps/int/1.1_genesKept.Rds")
	saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
	corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/drerio/nps/int/1.2_corrMat.Rds")
	saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}

if (args[[1]] == "amexicanus") {
  hypo.ast <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/AstMex_Hypo/AstMex_63k.rds")
	exprMat <- as.matrix(GetAssayData(hypo.ast))
	cellInfo <- data.frame(hypo.ast@meta.data)
	saveRDS(cellInfo, file="int/cellInfo.Rds")
	rm(hypo.ast)
	
	if (args[[2]] == "nps") {
	  org="mgi" # or hgnc, or dmel
	  dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/cisTarget_databases" # RcisTarget databases location
	  myDatasetTitle= paste("SCENIC analysis:", args[[1]], args[[2]], sep = " ") # choose a name for your analysis
	  data(defaultDbNames)
	  dbs <- defaultDbNames[[org]]
	  scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10)

	  scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
	  saveRDS(scenicOptions, file="int/scenicOptions.Rds")

	  counts <- rowSums(as.matrix(exprMat), na.rm = T)["hcrt"]
	  cells <- rowSums(as.matrix(exprMat) > 0, na.rm = T)["hcrt"]

	  genesKept <- geneFilteringv2(exprMat, scenicOptions=scenicOptions, minCountsPerGene=counts, minSamples=cells)
	  saveRDS(genesKept, file = "/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/amexicanus/nps/int/1.1_genesKept.Rds")

	  exprMat <- exprMat[genesKept, ]

	  tic("correlation")
	  runCorrelation(exprMat, scenicOptions)
	  toc()
	}
	
	genesKept <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/amexicanus/nps/int/1.1_genesKept.Rds")
	saveRDS(genesKept, file = "int/1.1_genesKept.Rds")
	corrMat <- readRDS("/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/amexicanus/nps/int/1.2_corrMat.Rds")
	saveRDS(corrMat, file = "int/1.2_corrMat.Rds")
}


## Initialize SCENIC settings (using mouse as stand in)

org="mgi" # or hgnc, or dmel
dbDir="/Volumes/BZ/Home/gizevo30/R_Projects/Cavefish_Paper/Seurat_v3_Integration/SCENIC/cisTarget_databases" # RcisTarget databases location
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
# Run the remaining steps using the wrapper functions:
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123


scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run

tic(paste("SCENIC step 1 with ", nrow(exprMat), " genes, ", ncol(exprMat), " cells, and ", length(TF_list), " TFs", sep = ""))
runSCENIC_1_coexNetwork2modules(scenicOptions)
toc()
