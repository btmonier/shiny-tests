  library(shiny)
  library(plotly)


  # User interface ----
  shinyApp(
    ui = fluidPage(
      tabsetPanel(
        tabPanel(
          "Welcome", icon = icon("bullhorn"),
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel = NULL,
            mainPanel = mainPanel(
              h1("Welcome to the ViDGER Pipeline"),
              h3("What in the world is this?"),
              p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer nec odio. Praesent libero. Sed cursus ante dapibus diam. Sed nisi. Nulla quis sem at nibh elementum imperdiet. Duis sagittis ipsum. Praesent mauris. Fusce nec tellus sed augue semper porta. Mauris massa. Vestibulum lacinia arcu eget nulla. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Curabitur sodales ligula in libero."),
              h3("Analyze"),
              p("Sed dignissim lacinia nunc. Curabitur tortor. Pellentesque nibh. Aenean quam. In scelerisque sem at dolor. Maecenas mattis. Sed convallis tristique sem. Proin ut ligula vel nunc egestas porttitor. Morbi lectus risus, iaculis vel, suscipit quis, luctus non, massa. Fusce ac turpis quis ligula lacinia aliquet. Mauris ipsum."),
              h3("Visualize"),
              p("Nulla metus metus, ullamcorper vel, tincidunt sed, euismod in, nibh. Quisque volutpat condimentum velit. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Nam nec ante. Sed lacinia, urna non tincidunt mattis, tortor neque adipiscing diam, a cursus ipsum ante quis turpis. Nulla facilisi. Ut fringilla. Suspendisse potenti. Nunc feugiat mi a tellus consequat imperdiet. Vestibulum sapien. Proin quam. Etiam ultrices. "),
              h3("Interact?"),
              p("Maybe...")
            )
          )
        ),
        tabPanel(
          "Analyze ", icon = icon("cogs"), 
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              selectInput(
                "Country", "Data", choices = "", selected = ""
              )
            ),
            mainPanel(
              htmlOutput("Attacks")
            )
          )
        ),
        tabPanel(
          "Visualize ", icon = icon("bar-chart"), 
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel(
              sliderInput(
                "year", "Year:", min = 1968, max = 2009, value = 2009, sep=''
              )
            ),
            mainPanel(
              fluidRow(
                column(7,  
                  plotlyOutput("")
                ),
                column(5, 
                  plotlyOutput("")
                )   
              )
            )
          )
        ),
        tabPanel(
          "Interact? ", icon = icon("binoculars"),
          fluid = TRUE,
          sidebarLayout(
            sidebarPanel = NULL,
            mainPanel = mainPanel(
              p("Interact with plotly here...")
            )
          )
        )
      )
    ),


  # Server logic ----
    server = function(input, output) {
      
    }
  )