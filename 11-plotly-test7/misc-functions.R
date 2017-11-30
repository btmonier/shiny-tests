#------------------------------------------------------------------------------
# Title:  ViDGER Shiny - Miscellaneous functions
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.24.17
#------------------------------------------------------------------------------

# Size factor extraction (for normalization)
getSizeFact <- function(rc.data) {
	geomMean <- function(x) { 
		prod(x)^(1 / length(x)) 
	}
	gm.mean <- apply(rc.data, 1, geomMean)
	gm.mean[gm.mean == 0] <- NA
	rc.data <- sweep(rc.data, 1, gm.mean, FUN = "/") 
	sf <- apply(rc.data, 2, median, na.rm = TRUE)
	return(sf)
}

# Get normalized counts from various object classes
getNormCounts <- function(rc.data) {
	if (class(rc.data) == "DESeqDataSet") {
		nc.data <- BiocGenerics::counts(rc.data, normalize = TRUE)
		return(nc.data)
	} else {
		sf <- getSizeFact(rc.data)
		nc.data <- sweep(rc.data, 2, sf, FUN = "/")
		return(nc.data)
	}
}

# Get normalized counts for selected gene (heatmap interactivity)
getGenes <- function(rc.data, id, coldata) {
	nc.data <- getNormCounts(rc.data)
	nc.data <- as.data.frame(nc.data[rownames(nc.data) == id, ])
	names(nc.data) <- "counts"
	dat.l <- list(coldata, nc.data)
	nc.data <- Reduce(
		merge, lapply(dat.l, function(x) data.frame(x, sample = row.names(x)))
	)
	return(nc.data)
}