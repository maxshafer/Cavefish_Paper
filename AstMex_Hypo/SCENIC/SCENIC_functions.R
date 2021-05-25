## Need to modify a couple of commands so that I can input my own custom TF list (or other list)

## Modify runGenie3 to use TF_list instead of the call to the Risc objects
runGenie3v2 <- function (exprMat, scenicOptions, nParts = 20, GO_list = GO_list, nTrees = 500, start = 1, weightMatrices = weightMatrices, ...) 
{
	require(tictoc)
    nCores <- getSettings(scenicOptions, "nCores")
    if (is.data.frame(exprMat)) {
        supportedClasses <- paste(gsub("AUCell_buildRankings,", 
            "", methods("AUCell_buildRankings")), collapse = ", ")
        supportedClasses <- gsub("-method", "", supportedClasses)
        stop("'exprMat' should be one of the following classes: ", 
            supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
    }
    if (any(table(rownames(exprMat)) > 1)) 
        stop("The rownames (gene id/name) in the expression matrix should be unique.")
    allTFs <- TF_list
    inputTFs <- allTFs[allTFs %in% rownames(exprMat)]
    percMatched <- length(inputTFs)/length(allTFs)
    if (getSettings(scenicOptions, "verbose")) 
        message("Using ", length(inputTFs), " TFs as potential regulators...")
    if (percMatched < 0.4) 
        warning("Only ", round(percMatched * 100), "% of the ", 
            length(allTFs), " TFs in the database were found in the dataset. Do they use the same gene IDs?\n")
    # genesSplit <- suppressWarnings(split(sort(rownames(exprMat)), 
        # 1:nParts))
    genesSplit <- suppressWarnings(split(sort(GO_list), 1:nParts))
    if (start == 1) {
    	weightMatrices <- list()
    } else {
    	weightMatrices <- weightMatrices
    }
    for (i in start:length(genesSplit)) {
        if (getSettings(scenicOptions, "verbose")) 
            message("Running GENIE3 part ", i)
        set.seed(getSettings(scenicOptions, "seed"))
        tic(paste("Running GENIE3 part ", i))
        weightMatrix <- GENIE3::GENIE3(exprMat, regulators = inputTFs, 
            nCores = nCores, targets = genesSplit[[i]], nTrees = nTrees)
        toc()
        fileName <- gsub(".Rds$", paste0("_part_", i, ".Rds"), 
            getIntName(scenicOptions, "genie3wm"))
        saveRDS(weightMatrix, file = fileName)
        weightMatrices[[i]] <- weightMatrix
    }
    linkList_list <- list()
    for (i in 1:length(genesSplit)) {
        weightMatrix <- weightMatrices[[i]]
        linkList_list[[i]] <- GENIE3::getLinkList(weightMatrix, 
            threshold = getSettings(scenicOptions, "modules/weightThreshold"))
    }
    rm(weightMatrices)
    linkList <- do.call(rbind, linkList_list)
    colnames(linkList) <- c("TF", "Target", "weight")
    linkList <- linkList[order(linkList[, "weight"], decreasing = TRUE), 
        ]
    saveRDS(linkList, file = getIntName(scenicOptions, "genie3ll"))
    if (getSettings(scenicOptions, "verbose")) 
        message("Finished running GENIE3.")
    invisible(linkList)
}

# Need to modify this function to remove the filtering based on the Risc Oject
# Change "genesKept <- genesLeft_minCells_inDatabases" to "genesKept <- genesLeft_minCells"

geneFilteringv2 <- function (exprMat, scenicOptions, minCountsPerGene = 3 * 0.01 * 
    ncol(exprMat), minSamples = ncol(exprMat) * 0.01) 
{
    outFile_genesKept <- NULL
    dbFilePath <- NULL
    if (class(scenicOptions) == "ScenicOptions") {
        dbFilePath <- getDatabases(scenicOptions)[[1]]
        outFile_genesKept <- getIntName(scenicOptions, "genesKept")
    }
    else {
        dbFilePath <- scenicOptions[["dbFilePath"]]
        outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
    }
    if (is.null(dbFilePath)) 
        stop("dbFilePath")
    if (is.data.frame(exprMat)) {
        supportedClasses <- paste(gsub("AUCell_buildRankings,", 
            "", methods("AUCell_buildRankings")), collapse = ", ")
        supportedClasses <- gsub("-method", "", supportedClasses)
        stop("'exprMat' should be one of the following classes: ", 
            supportedClasses, "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
    }
    if (any(table(rownames(exprMat)) > 1)) 
        stop("The rownames (gene id/name) in the expression matrix should be unique.")
    nCountsPerGene <- rowSums(exprMat, na.rm = T)
    nCellsPerGene <- rowSums(exprMat > 0, na.rm = T)
    message("Maximum value in the expression matrix: ", max(exprMat, 
        na.rm = T))
    message("Ratio of detected vs non-detected: ", signif(sum(exprMat > 
        0, na.rm = T)/sum(exprMat == 0, na.rm = T), 2))
    message("Number of counts (in the dataset units) per gene:")
    print(summary(nCountsPerGene))
    message("Number of cells in which each gene is detected:")
    print(summary(nCellsPerGene))
    message("\nNumber of genes left after applying the following filters (sequential):")
    genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > 
        minCountsPerGene)]
    message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", 
        minCountsPerGene)
    nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
    genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > 
        minSamples)]
    message("\t", length(genesLeft_minCells), "\tgenes detected in more than ", 
        minSamples, " cells")
    library(RcisTarget)
    motifRankings <- importRankings(dbFilePath)
    genesInDatabase <- colnames(getRanking(motifRankings))
    genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% 
        genesInDatabase)]
    message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
    genesKept <- genesLeft_minCells
    if (!is.null(outFile_genesKept)) {
        saveRDS(genesKept, file = outFile_genesKept)
        if (getSettings(scenicOptions, "verbose")) 
            message("Gene list saved in ", outFile_genesKept)
    }
    return(genesKept)
}



## Plotting functions

extractPlotData <- function(linkLists = dataset.lists, modules = dataset.modules, go.term = go.term, gene.name = gene.name) {
  
  linkLists <- lapply(linkLists, function(x) x[[go.term]][x[[go.term]]$Target == gene.name,])
  modules <- lapply(modules, function(x) x[[go.term]][x[[go.term]]$Target == gene.name, c(1,2,4)])
  
  for (i in 1:length(linkLists)) {
    rownames(linkLists[[i]]) <- linkLists[[i]]$TF
    linkLists[[i]]$corr <- modules[[i]]$corr[match(linkLists[[i]]$TF, modules[[i]]$TF)]
    linkLists[[i]]$corr[linkLists[[i]]$corr == 0] <- -1
    linkLists[[i]]$weight2 <- linkLists[[i]]$weight*linkLists[[i]]$corr
  }
  
  plot.data <- Reduce(function(x, y) merge(x, y, by = c("TF", "Target"), all = T), linkLists)
  colnames(plot.data) <- c("TF", "Target", unlist(lapply(names(linkLists), function(x) lapply(c("weight", "corr", "weight2"), function(y) c(paste(y, x, sep = "_"))))))
  plot.data[is.na(plot.data)] <- 0
  return(plot.data)
}

plotPlotData <- function(data = data, gene = gene, x = "surface", y = "cave", cutoff = NA, quantile = 0.95, log = F) {
  data <- data[[gene]]
  if (log) {
    data[[paste("weight", x, sep = "_")]] <- log(data[[paste("weight", x, sep = "_")]])
    data[[paste("weight", y, sep = "_")]] <- log(data[[paste("weight", y, sep = "_")]])
  }
  if (!is.na(cutoff)) {
    cutoff_x <- cutoff
    cutoff_y <- cutoff
  } else {
    cutoff_x <- quantile(data[[paste("weight", x, sep = "_")]][data[[paste("weight", x, sep = "_")]] != 0], probs = quantile, na.rm = T)
    cutoff_y <- quantile(data[[paste("weight", y, sep = "_")]][data[[paste("weight", y, sep = "_")]] != 0], probs = quantile, na.rm = T)
  }
  data$cutoff_x <- ifelse(data[[paste("weight", x, sep = "_")]] > cutoff_x, "yes", "no")
  data$cutoff_y <- ifelse(data[[paste("weight", y, sep = "_")]] > cutoff_y, "yes", "no")
  data$cutoff <- ifelse(data$cutoff_x == "yes" & data$cutoff_y == "yes", "both", ifelse(data$cutoff_x == "yes" & data$cutoff_y == "no", x, ifelse(data$cutoff_x == "no" & data$cutoff_y == "yes", y, "none")))
  
  data$cutoff <- factor(data$cutoff, levels = c(x, "both", "none", y))
  
  p <- ggplot(data, aes_string(x = paste("weight2", x, sep = "_"), y = paste("weight2", y, sep = "_"), colour = "cutoff")) + geom_point(size = 0.5) 
  p <- p + scale_colour_manual(values = c("#FF2800", "#FFB400", "black", "#06799F")) + theme_classic()
  p <- p + geom_text_repel(data = data[data$cutoff != "none",], aes_string(x = paste("weight2", x, sep = "_"), y = paste("weight2", y, sep = "_"), label = "TF"), size = 6, force = 5)
  p <- p + theme(axis.line = element_blank()) + geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
  if(log) {
    p <- p + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5)) + xlab(paste("log(Correlation with ", x, "-morph TFs)", sep = "")) + ylab(paste("log(Correlation with ", y, "-morph TFs)", sep = ""))
  } else {
    p <- p + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5)) + xlab(paste("Correlation with ", x, "-morph TFs", sep = "")) + ylab(paste("Correlation with ", y, "-morph TFs", sep = ""))
  }
  
  return(p)
}

calcDriftIndex <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	DI <- list()
	if (is.null(subset)) {
		DI <- lapply(seq_along(conserved), function(x) sqrt( abs( (1 - (length(conserved[[x]]) / length(species.1[[x]]))) * (1 - (length(conserved[[x]]) / length(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		DI <- lapply(seq_along(conserved), function(x) sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		DI <- lapply(seq_along(conserved), function(x) sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
	}
	}
	names(DI) <- names
	return(DI)
}

calcParalog <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, mart.1 = mart[[1]], mart.2 = mart[[2]], ngenes.1 = 32191, ngenes.2 = 25271, i = i) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	
	paralog.con <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(conserved[[i]], mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(conserved[[i]], mart.2$Gene.name)])
	paralog.1 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(setdiff(species.1[[i]], conserved[[i]]), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(setdiff(species.1[[i]], conserved[[i]]), mart.1$Gene.name)])
	paralog.2 <- union(mart.2$Zebrafish.paralogue.associated.gene.name[match(setdiff(species.2[[i]], conserved[[i]]), mart.2$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(setdiff(species.2[[i]], conserved[[i]]), mart.2$Gene.name)])
	paralog.union <- union(paralog.con, union(paralog.1, paralog.2))
	
	genes.1 <- setdiff(species.1[[i]], conserved[[i]])
	genes.2 <- setdiff(species.2[[i]], conserved[[i]])
	
	a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.2)])
	b1 <- length(genes.1) - a1
	c1 <- length(union(paralog.con, paralog.2)) - a1
	d1 <- ngenes.1 - b1

	a2 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.1)])
	b2 <- length(genes.1) - a2
	c2 <- length(union(paralog.con, paralog.1)) - a2
	d2 <- ngenes.2 - b2
	
	a3 <- length(conserved[[i]])
	b3 <- length(genes.1)
	c3 <- length(genes.2)
	d3 <- ngenes.1 - sum(b3, c3)
	
	vec <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2, a3=a3, b3=b3, c3=c3, d3=d3)
	return(vec)
}