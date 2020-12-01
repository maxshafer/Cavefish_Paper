# getBM_MS retrieves ensembl gene symbols given a chromosome and Start Stop positions
# Modification of biomaRt get_BM function
# Used to retrieve genes associated with FST windows

getBM_MS <- function(x, window = 25000) {
	require(biomaRt)
	CHROM <- x[1]
	START <- as.numeric(x[2])-window
	STOP <- as.numeric(x[3])+window
	genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name'), filters = c('chromosome_name', 'start', 'end'), values = list(CHROM, max(START, 0), STOP), mart = ensembl)
	if (nrow(genes) > 0)
		genes$SNP <- x[7]
	return(genes)
}

# Above is very slow, has to make new query every time
# getGenesEns takes genes (my_gene) from ensembl gene annotations (downloaded) (using the refGenome package)
getGenesEns <- function(x, gene_table = gtf2, window = 25000) {
	chrom <- as.character(x[1]) # change to as.character
	start <- as.numeric(x[2])-window
	stop <- as.numeric(x[3])+window
	cols <- c(1,3,4:5,9,22) # c(14, 17, 2, 4:5) # need to figure out which columns are needed from the old df
	gene_table <- gene_table[which(gene_table[1] == chrom), ] # 1 not 2
	gene_ids.1 <- gene_table[which(start > gene_table$start & start < gene_table$end), cols]
	gene_ids.2 <- gene_table[which(stop > gene_table$start & stop < gene_table$end), cols]
	gene_ids.3 <- gene_table[which(start < gene_table$start & stop > gene_table$end), cols]
	gene_ids <- rbind(gene_ids.1, gene_ids.2, gene_ids.3)
	if (nrow(gene_ids) > 0) {
		gene_ids$SNP <- as.numeric(x[7])
		gene_ids$CHROM <- as.character(x[1])
		gene_ids$BIN_START <- as.numeric(x[2])
		gene_ids$BIN_END <- as.numeric(x[3])
		gene_ids$N_VARIANTS <- as.numeric(x[4])
		gene_ids$MEAN_FST <- as.numeric(x[6])
		gene_ids$zscore <- as.numeric(x[10])
	}
	return(gene_ids)
}

# getGenesEns <- function(x, gene_table = my_gene, window = 25000) {
	# require(refGenome)
	# chrom <- x[1]
	# start <- as.numeric(x[2])-window
	# stop <- as.numeric(x[3])+window
	# fst <- seq(start, stop)
	# gene_ids <- apply(gene_table, 1, function(x) any(seq(x[4], x[5]) %in% fst)) # 30 mins per fst window?!
	# return(gene_ids)
# }

# readFST returns a data.frame of FST values
# Can use an fst cutoff value (0.1 for indels, 0.2 for snps), or z-transform fst values

readFST <- function(x, fisherZ = TRUE, variants = "windowed") {
	require(DescTools)
	fst.table <- read.table(x, header = T, stringsAsFactors = FALSE)
	fst.table <- fst.table[complete.cases(fst.table),]
	fst.table$SNP <- c(1:(nrow(fst.table)))
	# fst.table <- fst.table[fst.table$MEAN_FST > 0.15,]
	if (variants == "windowed") {
		fst.table$zscore <- (fst.table$MEAN_FST-mean(fst.table$MEAN_FST))/sd(fst.table$MEAN_FST)
	}
	if (variants == "perbase") {
		fst.table <- fst.table[fst.table$WEIR_AND_COCKERHAM_FST > 0.1, ]
		fst.table$zscore <- (fst.table$WEIR_AND_COCKERHAM_FST-mean(fst.table$WEIR_AND_COCKERHAM_FST))/sd(fst.table$WEIR_AND_COCKERHAM_FST)
	}
	if (fisherZ == TRUE) {
		fst.table$FisherZWeighted <- FisherZ(fst.table$WEIGHTED_FST)
		fst.table$FisherZMean <- FisherZ(fst.table$MEAN_FST)
		fst.table <- fst.table[!(fst.table$FisherZMean == Inf),]
		fst.table <- fst.table[!(fst.table$FisherZMean == -Inf),]
	}
	return(fst.table)
}

# Determine fdr cutoffs for weir.fst tables (returned by readFST)

fdrcutoff <- function(x, cutoff = 0.05) {
	require(fdrtool)
	require(dplyr)
	fdr <- fdrtool(x$zscore, plot = T, verbose = F)
	weir <- x[tibble(fdr[[3]]) < cutoff,]
	return(weir)
}