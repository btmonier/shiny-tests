#------------------------------------------------------------------------------
# Title:  Shiny Test 07 - Modularity - Tabs
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.26.17
#------------------------------------------------------------------------------

# Welcome page ----
tab.welcome <- tabPanel(
      title = "Welcome", icon = icon("hand-spock-o"),
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          HTML("<h1>Welcome to the ViDGER<sup>2</sup> Pipeline</h1>"),
          h3("What in the world is this?"),
          p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer nec odio. Praesent libero. Sed cursus ante dapibus diam. Sed nisi. Nulla quis sem at nibh elementum imperdiet. Duis sagittis ipsum. Praesent mauris. Fusce nec tellus sed augue semper porta. Mauris massa. Vestibulum lacinia arcu eget nulla. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Curabitur sodales ligula in libero."),
          h3("Analyze"),
          p("Sed dignissim lacinia nunc. Curabitur tortor. Pellentesque nibh. Aenean quam. In scelerisque sem at dolor. Maecenas mattis. Sed convallis tristique sem. Proin ut ligula vel nunc egestas porttitor. Morbi lectus risus, iaculis vel, suscipit quis, luctus non, massa. Fusce ac turpis quis ligula lacinia aliquet. Mauris ipsum."),
          h3("Visualize"),
          p("Nulla metus metus, ullamcorper vel, tincidunt sed, euismod in, nibh. Quisque volutpat condimentum velit. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Nam nec ante. Sed lacinia, urna non tincidunt mattis, tortor neque adipiscing diam, a cursus ipsum ante quis turpis. Nulla facilisi. Ut fringilla. Suspendisse potenti. Nunc feugiat mi a tellus consequat imperdiet. Vestibulum sapien. Proin quam. Etiam ultrices. "),
          h3("Interact?"),
          p("Yes.")
        )
      )
    )



# Submit and QC page ----
tab.submit <- tabPanel(
      title = "Submit and QC ", icon = icon("cogs"), 
      fluid = TRUE,
       sidebarLayout(
        sidebarPanel(
          h4("Parameters"),
          fileInput(
            inputId = "file1", 
            label = "Submit count data (CSV)",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            )
          ),
          fileInput(
            inputId = "file2", 
            label = "Submit sample information (CSV)",
            accept = c(
              "text/csv",
              "text/comma-separated-values,text/plain",
              ".csv"
            )
          ),
          textInput(
            "expdes",
            "Experimental design",
            value = ""
          ),
          textInput(
            "prefilt",
            "Filter cutoff",
            value = ""
          ),
          actionButton(inputId = "go", "Submit", icon = icon("magic")),
          uiOutput("opttran")
        ),
        mainPanel(
          br(),
          uiOutput("countsummary"),
          br(),
          # verbatimTextOutput("contents"),
          br(),
          uiOutput("countbox"),
          plotlyOutput("boxplot"),
          br(),
          uiOutput("counthist"),
          plotlyOutput("hist"),
          br(),
          uiOutput("counttotal"),
          plotlyOutput("barplot"),
          br(),
          uiOutput("countpca")
        )
      )
    )



# DEG Analysis ----
tab.deg <- tabPanel(
      title = "DEG Analysis ", icon = icon("bar-chart"), 
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel(
          h4("Parameters"),
          selectInput(
            inputId = "dge",
            label = "Choose analysis type",
            choices = c(
              "DESeq2" = "deseq",
              "edgeR" = "edger",
              "limma" = "limma"
            )
          ),
          uiOutput("degfact"),
          uiOutput("deglev1"),
          uiOutput("deglev2"),
          selectInput(
            inputId = "plottype",
            label = "Choose plot type",
            choices = c(
              "MA plot" = "maplot",
              "Volcano plot" = "volplot"
            )
          ),
          actionButton("godeg", "Submit", icon = icon("magic"))          
        ),
        mainPanel = mainPanel(
          br(),
          uiOutput("degcomp"),
          plotlyOutput("degplot"),
          br(),
          br(),
          DT::dataTableOutput("mytable"),

          uiOutput("downbutton"),
          downloadButton("downloadData", "Download Filtered Data"),
          br(),
          br(),
          br()
        )
      )
    )



# DEG ----
tab.help <- tabPanel(
      title = "Help ", icon = icon("question-circle"),
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          br(),
          p("Add a simple walkthrough?"),
          p("Add an FAQ?")
        )
      )
    )
