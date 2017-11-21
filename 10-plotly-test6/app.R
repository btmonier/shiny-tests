#------------------------------------------------------------------------------
# Title:  Shiny Test 08 - Crosstalk - Application
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.03.17
#------------------------------------------------------------------------------

# Package logic ----

## Set working directory (LOCAL ONLY)
# setwd("D:/Box Sync/misc-shiny-apps/10-plotly-test6")

## CRAN
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("crosstalk", "dplyr", "DT", "plotly", "shiny", "shinycssloaders",
			        "tidyr")
pack.man(packages)

## Bioconductor
source("https://bioconductor.org/biocLite.R")
if (!require("DESeq2")) biocLite("DESeq2")
if (!require("edgeR")) biocLite("edgeR")



# Sources ----
source("vidgerUI.R")
source("vidgerServer.R")



# Run shiny ----
shinyApp(vidgerUI, vidgerServer)