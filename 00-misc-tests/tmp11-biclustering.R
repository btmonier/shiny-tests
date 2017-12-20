#------------------------------------------------------------------------------
# Title:  Shiny Test 13 - Biclustering
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   12.14.17
#------------------------------------------------------------------------------

# Preamble

## Set working directory
setwd("D:/Box Sync/misc-shiny-apps/00-misc-tests")

## Load packages
library(QUBIC)
library(fields)

## Load data
cts <- read.csv("count-data.csv", header = TRUE, row.names = 1)
cts <- as.matrix(cts)
cts.log <- log2(cts + 1)


# Make biclusters

## Top variable genes
num <- 1000
cts.var <- order(rowVars(cts), decreasing = TRUE)
cts.var <- cts.log[cts.var, ]
cts.var <- cts.var[1:num, ]

## Run QUBIC algorithm
res <- biclust::biclust(cts.var, method = BCQU())

## Visualize
hmcols <- colorRampPalette(
	rev(
		c(
			"#D73027", "#FC8D59", "#FEE090", 
			"#FFFFBF", "#E0F3F8", "#91BFDB", 
			"#4575B4"
		)
	)
)(100)
par(mar = c(4, 5, 3, 5) + 0.1)
quheatmap(cts.var, res, number = c(1, 2), col = hmcols, showlabel = TRUE)