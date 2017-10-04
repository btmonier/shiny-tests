#-------------------------------------------------------------------------------
# Title:  Shiny - Vidger App Prototype
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   09.27.17
#-------------------------------------------------------------------------------

# Load packages ---
library(shiny)
library(vidger)
library(DESeq2)
library(edgeR)

# User interface ---
ui <- fluidPage(
  titlePanel("ViDGER Interactivity Test"),
  
  sidebarLayout(
    sidebarPanel(
      helpText("Let's make a box plot..."),
      
      selectInput("dataset", 
                  label = "Choose an example data set:", 
                  choices = c("df.cuff", "df.deseq", "df.edger"), 
                  selected = "df.cuff"),
      
      conditionalPanel(
        condition = "input.dataset == 'df.deseq'",
        textInput(inputId = "d.factor",
                  label = "DESeq2 Factor Level",
                  value = "cell")
      ),

      selectInput("ana.type", 
                  label = "Choose an analysis type:", 
                  choices = c("cuffdiff", "DESeq2", "edgeR"), 
                  selected = "cuffdiff"),

      radioButtons(inputId = "title", 
                   label = "Title?", 
                   choices = c("Yes", "No")),
      
      radioButtons(inputId = "legend", 
                   label = "Legend?", 
                   choices = c("Yes", "No")),
      
      radioButtons(inputId = "grid", 
                   label = "Grids?", 
                   choices = c("Yes", "No")),
      actionButton(inputId = "go", "Plot"),
      br(),
      br(),
      br(),
      br(),
      img(HTML('<center><img src="LabLogo.png" width="250"></center>'))
    ),
    mainPanel(plotOutput(outputId = "boxplot")))
)


# Server logic ---
server <- function(input, output) {
  output$boxplot <- renderPlot({
    
    input$go
    
    data <- isolate(switch(input$dataset,
                   "df.cuff" = df.cuff,
                   "df.deseq" = df.deseq,
                   "df.edger" = df.edger))

    isolate(
      if (input$d.factor == "") {
        d.factor <- NULL
      } else {
        d.factor <- input$d.factor
      }
    )
    
    type <- isolate(switch(input$ana.type,
                   "cuffdiff" = "cuffdiff",
                   "DESeq2" = "deseq",
                   "edgeR" = "edger"))

    title <- isolate(switch(input$title,
                    "Yes" = TRUE,
                    "No" = FALSE))

    legend <- isolate(switch(input$legend,
                     "Yes" = TRUE,
                     "No" = FALSE))

    grid <- isolate(switch(input$grid,
                   "Yes" = TRUE,
                   "No" = FALSE))

    vsBoxPlot(data, d.factor, type, title, legend, grid)

  })
}


# Run app ----
shinyApp(ui, server)