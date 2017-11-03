#------------------------------------------------------------------------------
# Title:  Shiny Test 07 - Modularity - Application
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.19.17
#------------------------------------------------------------------------------

# Package logic ----

## CRAN
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("crosstalk", "shiny", "DT", "tidyr", "plotly")
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