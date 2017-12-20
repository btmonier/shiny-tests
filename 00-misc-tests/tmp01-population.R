#------------------------------------------------------------------------------
# Title:  selectInput population test
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   09.29.17
#------------------------------------------------------------------------------

library("shiny")
ui <- fluidPage(
  titlePanel("selectInput population test"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        'file1', 
        'Choose CSV file',
        accept = c('.csv','.tsv')
      ),
      textInput(
        "text",
        "Define factor",
        ""
      ),
       
      # tags$script(
      #   '$( "#file1" ).on( "click", function() {
      #   this.value = null; 
      #   });'
      # ),
      uiOutput("opt.x")
    ),
    mainPanel(
      tableOutput("your_data")
    )
  )
)
server <- function(input, output) {
  
  tableData <- reactive({
    inFile <- input$file1
    if (!is.null(inFile)){
      read.csv(inFile$datapath)
    } else {
      NULL
    }
  })

  output$opt.x <- renderUI({
    tmp <- tableData()
    validate(
      need(input$text == names(tmp), "Please enter correct factor")
    )
        
    selectInput("xcolumn", "Your factor levels", tmp[, input$text])
  })
  
  output$your_data <- renderTable({
    head(tableData())
  })
}
shinyApp(ui = ui, server = server)