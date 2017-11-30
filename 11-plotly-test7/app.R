#------------------------------------------------------------------------------
# Title:  Shiny Test 11 - Crosstalk - Application (With heatmaps)
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.24.17
#------------------------------------------------------------------------------

# Package logic ----

## Set working directory (FOR LOCAL TESTING ONLY)
# setwd("D:/Box Sync/misc-shiny-apps/11-plotly-test7")

## CRAN
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("crosstalk", "dplyr", "DT", "plotly", "shiny", "shinyBS", 
		          "shinycssloaders", "tibble", "tidyr")
pack.man(packages)

## Bioconductor
source("https://bioconductor.org/biocLite.R")
if (!require("DESeq2")) biocLite("DESeq2")
if (!require("edgeR")) biocLite("edgeR")
if (!require("limma")) biocLite("limma")

## Example data files
f1 <- as.matrix(read.csv("count-data.csv", header = TRUE, row.names = 1))
f2 <- read.csv("col-data.csv", header = TRUE, row.names = 1)



# Sources ----
source("vidgerUI.R")
source("vidgerServer.R")



# Run shiny ----
shinyApp(vidgerUI, vidgerServer)