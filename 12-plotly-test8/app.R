#------------------------------------------------------------------------------
# Title:  Shiny Test 12 - Experimental Designs - Application
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   12.08.17
#------------------------------------------------------------------------------

# Package logic ----

## Set working directory (FOR LOCAL TESTING ONLY)
# setwd("D:/Box Sync/misc-shiny-apps/12-plotly-test8")

## CRAN
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("crosstalk", "dplyr", "DT", "gtools", "plotly", "shiny", 
							"shinyBS", "shinycssloaders", "shinythemes", "tibble", "tidyr",
							"Rcpp", "Hmisc", "ggplot2", "geneplotter", "locfit")
pack.man(packages)

## Bioconductor
source("https://bioconductor.org/biocLite.R")
if (!require("DESeq2")) biocLite("DESeq2")
if (!require("edgeR")) biocLite("edgeR")
if (!require("limma")) biocLite("limma")
if (!require("QUBIC")) biocLite("QUBIC")



# Sources ----
source("vidgerUI.R")
source("vidgerServer.R")



# Run shiny ----
shinyApp(vidgerUI, vidgerServer)