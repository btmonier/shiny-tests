#------------------------------------------------------------------------------
# Title:  ViDGER DESeq Hack Functions
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.17.17
#------------------------------------------------------------------------------


# Boxplot data extraction
getDeseqBox <- function(data, d.factor, tran) {

	if(is.null(d.factor)) {
	  stop(
	  	"This appears to be a DESeq object. Please state d.factor variable."
	  )
	}
	if (tran == "rlog") {
	  tmp <- rlog(dds, blind = FALSE)
	  tmp <- assay(tmp)
	  lab <- "rlog(counts)"
	} else if (tran == "vst") {
	  tmp <- vst(dds, blind = FALSE)
	  tmp <- assay(tmp)
	  lab <- "vst(counts)"
	} else if (tran == "counts") {
	  tmp <- assay(dds)
	  lab <- "raw counts"
	} else if (tran == "log") {
	  tmp <- normTransform(dds)
	  tmp <- assay(tmp)
	  lab <- "log2(counts + 1)"
	}
	dat1 <- as.data.frame(colData(dds))
	dat2 <- as.data.frame(tmp)
	nam <- as.vector(unique(dat1[[d.factor]]))
	ls.nam <- list()
	ls.mean <- list()
	for (i in nam) {
	  ls.nam[[i]] <- row.names(dat1[which(dat1[d.factor] == i), ])
	  for (j in 1:length(ls.nam)){
	    ls.mean[[j]] <- rowMeans(dat2[, ls.nam[[j]]])
	  }
	}
	names(ls.mean) <- sapply(nam, paste)
	dat3 <- as.data.frame(ls.mean)
	dat3 <- tidyr::gather(as.data.frame(dat3))
	dat3$key <- as.factor(dat3$key)
	return(dat3)
}


# PlotCounts data extraction
getDeseqGeneCount <- function(data, d.factor, tran, gene) {
	if(is.null(d.factor)) {
	  stop(
	  	"This appears to be a DESeq object. Please state d.factor variable."
	  )
	}
	if (tran == "rlog") {
	  tmp <- rlog(dds, blind = FALSE)
	  tmp <- assay(tmp)
	  lab <- "rlog(counts)"
	} else if (tran == "vst") {
	  tmp <- vst(dds, blind = FALSE)
	  tmp <- assay(tmp)
	  lab <- "vst(counts)"
	} else if (tran == "counts") {
	  tmp <- assay(dds)
	  lab <- "raw counts"
	} else if (tran == "log") {
	  tmp <- normTransform(dds)
	  tmp <- assay(tmp)
	  lab <- "log2(counts + 1)"
	}
	# tmp <- tmp[which(rownames(tmp) == gene), ]
	# tmp <- tidyr::gather(as.data.frame(t(tmp)))

	dat1 <- as.data.frame(colData(dds))
	dat2 <- as.data.frame(tmp)
	nam <- as.vector(unique(dat1[[d.factor]]))
	ls.nam <- list()
	ls.mean <- list()
	for (i in nam) {
	  ls.nam[[i]] <- row.names(dat1[which(dat1[d.factor] == i), ])
	  for (j in 1:length(ls.nam)){
	    ls.mean[[j]] <- dat2[, ls.nam[[j]]]
	  }
	}
	names(ls.mean) <- sapply(nam, paste)
	dat3 <- as.data.frame(ls.mean)
	dat3 <- dat3[which(rownames(dat3) == gene), ]
	dat3 <- tidyr::gather(as.data.frame(dat3))
	dat3 <- tidyr::separate(dat3, key, into = c(d.factor, "sample"), 
							sep = "\\.(?=[^\\.]*$)")
	dat3[[1]] <- as.factor(dat3[[1]])
	dat3[[2]] <- as.factor(dat3[[2]])
	return(dat3)

}


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