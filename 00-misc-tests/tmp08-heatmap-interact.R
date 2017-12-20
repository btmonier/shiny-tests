#------------------------------------------------------------------------------
# Title:  Heatmap interactions
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.22.17
#------------------------------------------------------------------------------

# Load packages
library(plotly)
library(shiny)
library(tibble)

setwd("D:/Box Sync/misc-shiny-apps/00-misc")

# compute a correlation matrix (run tmp3 first!)

## Create array (non DESeq S4 object)
dds <- pas.dds
res <- results(pas.dds, tidy = TRUE)
res <- as_tibble(res)
# vsd <- vst(pas.dds, blind = FALSE)
# rld <- rlog(pas.dds, blind = FALSE)
ntd <- normTransform(dds)

## Heat map (plotly)
num <- 25
heat <- assay(ntd)[arrange(res, padj, pvalue)$row[1:num], ]
heat <- t(heat)
heat <- scale(heat)
heat <- t(heat)

ui <- fluidPage(
  mainPanel(
    plotlyOutput("heat"),
    plotlyOutput("scatterplot")
  ),
  verbatimTextOutput("selection")
)

server <- function(input, output, session) {
  output$heat <- renderPlotly({
    plot_ly(
      x = 1:num,
      y = rownames(heat),
      z = heat,
      type = "heatmap",
      source = "heatplot"
    ) %>%
    layout(
      xaxis = list(title = ""),
      yaxis = list(title = "")
    )
  })
  
  output$selection <- renderPrint({
    s <- event_data("plotly_click", source = "heatplot")
    if (length(s) == 0) {
      "Click on a cell in the heatmap to display a scatterplot"
    } else {
      cat("You selected: \n\n")
      as.list(s)
    }
  })

  output$scatterplot <- renderPlotly({
    s <- event_data("plotly_click", source = "heatplot")
    source("tmp5-vidger-functions.R")
    test <- getGenes(
      rc.data = pas.dds, 
      id = s[["y"]],
      coldata = pas.coldata
    )

    plot_ly(
      data = test,
      type = "scatter",
      mode = "markers",
      x = ~condition,
      y = ~counts,
      color = ~condition
    )    
  })
}

shinyApp(ui, server, options = list(display.mode = "showcase"))
