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

extractPlotData <- function(surface = GENIE3_linkList.surface, cave = GENIE3_linkList.cave, cor.surface = tfModules_asDF.surface, cor.cave = tfModules_asDF.cave, np = np) {
  surface <- surface[surface$Target == np,]
  cave <- cave[cave$Target == np,]
  
  cor.cave <- unique(cor.cave[cor.cave$Target == np, c(1,2,4)])
  cor.surface <- unique(cor.surface[cor.surface$Target == np, c(1,2,4)])
  
  surface$corr <- cor.surface$corr[match(surface$TF, cor.surface$TF)]
  cave$corr <- cor.cave$corr[match(cave$TF, cor.cave$TF)]
  # intersect <- intersect(row.names(surface), row.names(cave))
  # surface <- surface[intersect,]
  # cave <- cave[intersect,]
  # plot.data <- cbind(surface, cave)
  
  plot.data <- merge(cave, surface, by = "TF", all = T)
  plot.data[,c(3,6)][is.na(plot.data[,c(3,6)])] <- 0
  colnames(plot.data) <- c("TF", "Target_cave", "weight_cave", "corr_cave", "Target_surface", "weight_surface", "corr_surface")
  return(plot.data)
}

plotPlotData <- function(data = data, gene = gene, cutoff = NA, quantile = 0.95, log = F) {
	data <- data[[gene]]
	if (log) {
		data$weight_cave <- log(data$weight_cave)
		data$weight_surface <- log(data$weight_surface)
	}
	if (!is.na(cutoff)) {
		cutoff_surface <- cutoff
		cutoff_cave <- cutoff
	} else {
		cutoff_surface <- quantile(data$weight_surface[data$weight_surface != 0], probs = quantile, na.rm = T)
		cutoff_cave <- quantile(data$weight_cave[data$weight_cave != 0], probs = quantile, na.rm = T)
	}
	data$cutoff_surface <- ifelse(data$weight_surface > cutoff_surface, "yes", "no")
	data$cutoff_cave <- ifelse(data$weight_cave > cutoff_cave, "yes", "no")
	data$cutoff <- ifelse(data$cutoff_surface == "yes" & data$cutoff_cave == "yes", "both", ifelse(data$cutoff_surface == "yes" & data$cutoff_cave == "no", "surface", ifelse(data$cutoff_surface == "no" & data$cutoff_cave == "yes", "cave", "none")))
	data$cutoff <- factor(data$cutoff, levels = c("both", "cave", "none", "surface"))
	p <- ggplot(data, aes(x = weight_cave, y = weight_surface, colour = cutoff)) + geom_point(size = 0.1) 
	p <- p + scale_colour_manual(values = c("#440154FF", "#FDE725FF", "black", "#21908CFF")) + theme_classic()
	p <- p + geom_text_repel(data = data[data$cutoff != "none",], aes(x = weight_cave, y = weight_surface, label = TF), size = 8/2.856, force = 5)
	if(log) {
		p <- p + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 10)) + xlab("log(Random forest weight \n Mexican tetra cave-morph TFs)") + ylab("log(Random forest weight \n Mexican tetra surface-morph TFs)")
	} else {
		p <- p + ggtitle(gene) + theme(plot.title = element_text(hjust = 0.5, face = "italic", size = 10)) + xlab("Random forest weight \n Mexican tetra cave-morph TFs") + ylab("Random forest weight \n Mexican tetra surface-morph TFs")
	}
	
	return(p)
}

calcSimilarityIndex <- function(conserved = conserved, species.1 = species.1, species.2 = species.2, subset = NULL, invert = FALSE) {
	names <- Reduce(intersect, list(names(conserved), names(species.1), names(species.2)))
	conserved <- conserved[names]
	species.1 <- species.1[names]
	species.2 <- species.2[names]
	DI <- list()
	if (is.null(subset)) {
		DI <- lapply(seq_along(conserved), function(x) 1- sqrt( abs( (1 - (length(conserved[[x]]) / length(species.1[[x]]))) * (1 - (length(conserved[[x]]) / length(species.2[[x]]))) ) ))
	} else {
		if (invert == FALSE) {
		DI <- lapply(seq_along(conserved), function(x) 1- sqrt( abs( (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset]) / length(row.names(species.1[[x]])[row.names(species.1[[x]]) %in% subset]))) * (1 - (length(row.names(conserved[[x]])[row.names(conserved[[x]]) %in% subset])) / length(row.names(species.2[[x]])[row.names(species.2[[x]]) %in% subset])) ) ))
	} else {
		DI <- lapply(seq_along(conserved), function(x) 1 -sqrt( abs( (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)]) / length(row.names(species.1[[x]])[!(row.names(species.1[[x]]) %in% subset)]))) * (1 - (length(row.names(conserved[[x]])[!(row.names(conserved[[x]]) %in% subset)])) / length(row.names(species.2[[x]])[!(row.names(species.2[[x]]) %in% subset)])) ) ))
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
	paralog.2 <- union(mart.1$Zebrafish.paralogue.associated.gene.name[match(setdiff(species.2[[i]], conserved[[i]]), mart.1$Gene.name)], mart.2$Cave.fish.paralogue.associated.gene.name[match(setdiff(species.2[[i]], conserved[[i]]), mart.2$Gene.name)])
	paralog.union <- union(paralog.con, union(paralog.1, paralog.2))
	
	paralog.1 <- paralog.1[paralog.1 %in% mart.2$Gene.name]
	paralog.2 <- paralog.2[paralog.2 %in% mart.1$Gene.name]
	
	genes.1 <- setdiff(species.1[[i]], conserved[[i]])
	genes.2 <- setdiff(species.2[[i]], conserved[[i]])
	
	a1 <- length(genes.1[genes.1 %in% union(paralog.con, paralog.2)])
	b1 <- length(genes.1) - a1
	c1 <- length(union(paralog.con, paralog.2)) - a1
	d1 <- ngenes.1 - b1

	a2 <- length(genes.2[genes.2 %in% union(paralog.con, paralog.1)])
	b2 <- length(genes.2) - a2
	c2 <- length(union(paralog.con, paralog.1)) - a2
	d2 <- ngenes.2 - b2
	
	a3 <- length(conserved[[i]])
	b3 <- length(genes.1)
	c3 <- length(genes.2)
	d3 <- ngenes.1 - sum(b3, c3)
	
	vec <- c(a1=a1, c1=c1, b1=b1, d1=d1, a2=a2, c2=c2, b2=b2, d2=d2, a3=a3, b3=b3, c3=c3, d3=d3)
	return(vec)
}