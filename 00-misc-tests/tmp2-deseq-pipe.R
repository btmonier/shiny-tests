#------------------------------------------------------------------------------
# Title:  Shiny - DESeq2 Workflow
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   09.28.17
#------------------------------------------------------------------------------

# Load packages ----
library(shiny)
library(DESeq2)


# User interface ----
ui <- fluidPage(
  titlePanel("DESeq2 Pipeline"),
  
  sidebarLayout(
    sidebarPanel(
      fileInput("file1", "Submit gene count data (CSV)",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      fileInput("file2", "Submit annotation data (CSV)",
                accept = c(
                  "text/csv",
                  "text/comma-separated-values,text/plain",
                  ".csv")
      ),
      textInput(
        "expdes",
        "Experimental design:",
        value = ""
      ),
      radioButtons("filt1", "Pre-filter data?", choices = c("Yes", "No")),
      actionButton("go","Submit")
      # br(),
      # br(),
      # br(),
      # br()
      # img(HTML('<center><img src="LabLogo.png" width="250"></center>'))
    ),
    mainPanel(
      h4("DESeq2 results summary"),
      verbatimTextOutput("contents"),
      br(),
      h4("DESeq2 plot results (MA)"),
      plotOutput("plot")
    )
  )
)

# Server logic ----
server <- function(input, output) {
  ddsout <- reactive({
    if (input$go == 0) {
      return()
    } else {
      isolate({
        cts <- input$file1
        coldata <- input$file2
        cts <- as.matrix(read.csv(cts$datapath, header = TRUE, row.names = 1))
        coldata <- read.csv(coldata$datapath, header = TRUE, row.names = 1)
        cts <- cts[, rownames(coldata)]
        
        if (input$expdes == "") {
          expdes <- NULL
        } else {
          expdes <- input$expdes
        }
        
        dds <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = as.formula(expdes))
        dds <- DESeq(dds)
        
        reslfc <- lfcShrink(dds = dds, coef = 2, res = results(dds))
        
        return(list(dds, reslfc))  
      })
        
    }
  })
  
  output$contents <- renderPrint({
    if (input$go == 0) {
      return()
    } else {
      isolate({
        if (input$filt1 == "Yes") {
          dds <- ddsout()[[1]]
          dds <- dds[ rowSums(counts(dds)) > 1, ]
        } else {
          dds <- ddsout()[[1]]
          summary(dds)
        }  
      })
      
    }
    summary(results(dds))
  })  
  
  output$plot <- renderPlot({
    if (input$go == 0) {
      return()
    } else {
      isolate({
        reslfc <- ddsout()[[2]]
        plotMA(reslfc,  ylim=c(-2,2))
      })  
    }
  })  
  
}

shinyApp(ui, server)

