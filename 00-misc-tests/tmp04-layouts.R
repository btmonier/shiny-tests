#------------------------------------------------------------------------------
# Title:  DESeq2 Console Backchannel
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.06.17
#------------------------------------------------------------------------------

# Load packages ----
library(shiny)
library(plotly


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
            p("Lorem ipsum dolor sit amet, consectetur adipiscing elitInteger nec odio. Praesent libero. Sed cursus ante dapibudiam. Sed nisi. Nulla quis sem at nibh elementum imperdietDuis sagittis ipsum. Praesent mauris. Fusce nec tellus seaugue semper porta. Mauris massa. Vestibulum lacinia arcu egenulla. Class aptent taciti sociosqu ad litora torquent peconubia nostra, per inceptos himenaeos. Curabitur sodaleligula in libero."),
            h3("Analyze"),
            p("Sed dignissim lacinia nunc. Curabitur tortor. Pellentesqunibh. Aenean quam. In scelerisque sem at dolor. Maecenamattis. Sed convallis tristique sem. Proin ut ligula vel nunegestas porttitor. Morbi lectus risus, iaculis vel, suscipiquis, luctus non, massa. Fusce ac turpis quis ligula lacinialiquet. Mauris ipsum."),
            h3("Visualize"),
            p("Nulla metus metus, ullamcorper vel, tincidunt sed, euismoin, nibh. Quisque volutpat condimentum velit. Class aptentaciti sociosqu ad litora torquent per conubia nostra, peinceptos himenaeos. Nam nec ante. Sed lacinia, urna notincidunt mattis, tortor neque adipiscing diam, a cursus ipsuante quis turpis. Nulla facilisi. Ut fringilla. Suspendisspotenti. Nunc feugiat mi a tellus consequat imperdietVestibulum sapien. Proin quam. Etiam ultrices. "),
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
  )


# Server logic ----
  server = function(input, output) {
    
  }
)