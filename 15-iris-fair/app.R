#---------------------------------------------------------------------
# Title:         IRIS - Shiny Application
# Author:        Brandon Monier
# Created:       2018-01-26 11:29:39 CDT
# Last Modified: 2018-05-09 13:31:20 CDT
#---------------------------------------------------------------------

# Package logic ----

## Set working directory (FOR LOCAL TESTING ONLY)
# setwd("D:/Box Sync/misc-github-repos/shiny-tests/15-iris-fair")

## CRAN
source("iris-functions.R")
packages <- c(
	"crosstalk", "dplyr", "DT", "gtools", "plotly", "shiny", "plyr",
	"shinyBS", "shinycssloaders", "shinythemes", "tibble", "tidyr",
	"Rcpp", "Hmisc", "ggplot2", "geneplotter", "locfit", "GGally", 
	"pheatmap",	"reshape2", "backports", "digest", "fields", "psych",
	"stringr", "tools", "openxlsx"
)
pack.man(packages)

## Bioconductor
source("https://bioconductor.org/biocLite.R")
if (!require("DESeq2")) biocLite("DESeq2")
if (!require("edgeR")) biocLite("edgeR")
if (!require("limma")) biocLite("limma")
if (!require("QUBIC")) biocLite("QUBIC")



# Sources ----
source("irisUI.R")
source("irisServer.R")
source("iris-xlsx.R")



# Run shiny ----
shinyApp(irisUI, irisServer)
