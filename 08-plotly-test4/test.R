library(shiny)

ui <- fluidPage(
  titlePanel("selectInput population test"),
  sidebarLayout(
    sidebarPanel(
      fileInput(
        'file1', 
        'Choose CSV file',
        accept = c('.csv','.tsv')
      ),
      # uiOutput("opt.x"),
      uiOutput("opt.y"),
      uiOutput("opt.z"),
      actionButton("godge", "Submit", icon = icon("magic"))
    ),
    mainPanel(
      tableOutput("your_data"),
      br(),
      uiOutput("header_panel"),
      plotOutput("scatter")
    )
  )
)
server <- function(input, output) {
  
  tableData <- reactive({
    inFile <- input$file1
    if (!is.null(inFile)) {
      read.csv(inFile$datapath)
    } else {
      NULL
    }
  })

  # output$opt.x <- renderUI({
  #   tmp <- tableData()
  #   selectInput("factor", "Define factor", colnames(tmp))

  # })

  output$opt.y <- renderUI({
    tmp <- tableData()
    selectInput("lev1", "Define level 1", colnames(tmp))
  })

  output$opt.z <- renderUI({
    tmp <- tableData()
    selectInput("lev2", "Define level 2", colnames(tmp))
  })
  
  output$your_data <- renderTable({
    head(tableData())
  })

  output$header_panel <- renderUI({
    if(input$godge == 0) {
      return()
    } else {
      h4("A scatter plot...")
    }
  })

  scatterplot <- eventReactive(input$godge, {
    tmp <- tableData()
    x <- tmp[[input$lev1]]
    y <- tmp[[input$lev2]]
    plot(x, y)
  })

  output$scatter <- renderPlot({
    scatterplot()
  })


  # output$scatter <- renderPlot(eventReactive(input$godge, {
    
  #     tmp <- tableData()

  #     plot(x, y)
  #   }))

    # validate(
    #   need(input$godge != 0, "")
    # )
    # tmp <- tableData()
    # plot(x, y)
  # })

}
shinyApp(ui = ui, server = server)