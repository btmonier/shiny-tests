#------------------------------------------------------------------------------
# Title:  Crosstalk Experiments
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.27.17
#------------------------------------------------------------------------------

# Load packages
library(crosstalk)
library(d3heatmap)
library(d3scatter)
library(leaflet)
library(plotly)
library(shiny)
library(DT)



# Load data
# setwd("D:/Box Sync/misc-shiny-apps/00-misc-tests")
datobj <- read.csv("deseq-test.csv", header = TRUE, row.names = 1)
datobj <- datobj[, -c(6, 7)]
# datobj <- round(datobj, 3)
# datobj$padj[is.na(datobj$padj)] <- 1


# Shiny server crosstalk 2
ui <- fluidPage(
  h2("Plotly & DT Test 2", align = "center"),
  selectInput(
    "plottype", 
    "Choose plot type:", 
    choices = c("MA plot" = "maplot", "Volcano plot" = "voplot"),
    selected = "voplot"
  ),
  textInput("padj", "Adjusted p-value less than:", value = "1"),
  textInput("lfc", "Absolute log2 foldchange greater than:", value = "0"),
  plotlyOutput("scatterplot"),
  br(),
  DT::dataTableOutput("table"),
  br()
)

server <- function(input, output) {
  
  datout <- reactive({
    datobj <- datobj[which(datobj$padj <= input$padj), ]
    datobj <- datobj[which(abs(datobj$log2FoldChange) >= input$lfc), ]
    datobj <- datobj %>% tibble::rownames_to_column()
  })

  share <- reactive({
    SharedData$new(datout())
  })

  output$scatterplot <- renderPlotly({
    s <- input$table_rows_selected

    if (input$plottype == "maplot") {
      if (!length(s)) {
        p <- share() %>%
          plot_ly(x = ~log10(baseMean), y = ~log2FoldChange, mode = "markers", color = I("darkgray"), name = "Unfiltered") %>%
          layout(showlegend = TRUE) %>% 
          highlight("plotly_selected", color = I("royalblue1"), selected = attrs_selected(name = "Filtered"))
      } else if (length(s)) {
        pp <- datout() %>%
          plot_ly() %>% 
          add_trace(x = ~log10(baseMean), y = ~log2FoldChange, mode = "markers", color = I("darkgray"), name = "Unfiltered") %>%
          layout(showlegend = TRUE)

        # selected data
        pp <- add_trace(pp, data = datout()[s, , drop = FALSE], x = ~log10(baseMean), y = ~log2FoldChange, mode = "markers",
                        color = I("royalblue1"), size = I(10), name = "Filtered")
      }
    } else if (input$plottype == "voplot") {
      if (!length(s)) {
        p <- share() %>%
          plot_ly(x = ~log2FoldChange, y = ~-log10(pvalue), mode = "markers", color = I("darkgray"), name = "Unfiltered") %>%
          layout(showlegend = TRUE) %>% 
          highlight("plotly_selected", color = I("royalblue1"), selected = attrs_selected(name = "Filtered"))
      } else if (length(s)) {
        pp <- datout() %>%
          plot_ly() %>% 
          add_trace(x = ~log2FoldChange, y = ~-log10(pvalue), mode = "markers", color = I("darkgray"), name = "Unfiltered") %>%
          layout(showlegend = TRUE)

        # selected data
        pp <- add_trace(pp, data = datout()[s, , drop = FALSE], x = ~log2FoldChange, y = ~-log10(pvalue), mode = "markers",
                        color = I("royalblue1"), size = I(10), name = "Filtered")
      }      
    }
  })

  # highlight selected rows in the table
  output$table <- DT::renderDataTable({
    datout()
    m2 <- datout()[share()$selection(),]
    dt <- DT::datatable(datout())
    if (NROW(m2) == 0) {
      dt
    } else {
      DT::formatStyle(dt, "rowname", target = "row",
                      color = DT::styleEqual(m2$rowname, rep("white", length(m2$rowname))),
                      backgroundColor = DT::styleEqual(m2$rowname, rep("darkgray", length(m2$rowname))))
    }
  })
}

shinyApp(ui, server)