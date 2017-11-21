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
      actionButton(inputId = "go", "Submit"),
      br(),
      br(),
      br(),
      br(),
      img(HTML('<center><img src="LabLogo.png" width="250"></center>'))
    ),
    mainPanel(
      h4("DESeq2 results summary"),
      verbatimTextOutput("contents")
    )
  )
)

server <- function(input, output) {
  output$contents <- renderPrint({
    input$go
    cts <- isolate(input$file1)
    coldata <- isolate(input$file2)
    
    if (is.null(cts))
      return(NULL)
    
    cts <- as.matrix(read.csv(cts$datapath, header = TRUE, row.names = 1))
    coldata <- read.csv(coldata$datapath, header = TRUE, row.names = 1)
    
    cts <- cts[, rownames(coldata)]
    
    isolate(
      if (input$expdes == "") {
        expdes <- NULL
      } else {
        expdes <- input$expdes
      }
    )

    dds <- DESeqDataSetFromMatrix(countData = cts,
                                  colData = coldata,
                                  design = as.formula(expdes))
    
    isolate(
      if (input$filt1 == "Yes") {
        dds <- dds[ rowSums(counts(dds)) > 1, ]
      } else {
        dds <- dds
      }
    )
    
    dds <- DESeq(dds)
    res <- results(dds)
    head(res)
    
  })
  

}

shinyApp(ui, server)

