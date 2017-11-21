#------------------------------------------------------------------------------
# Title:  Shiny Test 08 - Crosstalk - Tabs
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.03.17
#------------------------------------------------------------------------------

# Welcome page ----
tab.welcome <- tabPanel(
  title = "Welcome", icon = icon("hand-spock-o"),
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      HTML("<h1>Welcome to the ViDGER<sup>2</sup> Pipeline</h1>"),
      h3("About"),
      p("Lorem ipsum dolor sit amet, consectetur adipiscing elit. Integer nec odio. Praesent libero. Sed cursus ante dapibus diam. Sed nisi. Nulla quis sem at nibh elementum imperdiet. Duis sagittis ipsum. Praesent mauris. Fusce nec tellus sed augue semper porta. Mauris massa. Vestibulum lacinia arcu eget nulla. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Curabitur sodales ligula in libero."),
      h3("Submit and QC"),
      p("Sed dignissim lacinia nunc. Curabitur tortor. Pellentesque nibh. Aenean quam. In scelerisque sem at dolor. Maecenas mattis. Sed convallis tristique sem. Proin ut ligula vel nunc egestas porttitor. Morbi lectus risus, iaculis vel, suscipit quis, luctus non, massa. Fusce ac turpis quis ligula lacinia aliquet. Mauris ipsum."),
      h3("DEG Analysis"),
      p("Nulla metus metus, ullamcorper vel, tincidunt sed, euismod in, nibh. Quisque volutpat condimentum velit. Class aptent taciti sociosqu ad litora torquent per conubia nostra, per inceptos himenaeos. Nam nec ante. Sed lacinia, urna non tincidunt mattis, tortor neque adipiscing diam, a cursus ipsum ante quis turpis. Nulla facilisi. Ut fringilla. Suspendisse potenti. Nunc feugiat mi a tellus consequat imperdiet. Vestibulum sapien. Proin quam. Etiam ultrices. "),
      h3("Help"),
      p("Need a walkthrough? Look through this tab for detailed information about the pipeline.")
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
      actionButton(inputId = "goqc", "Submit", icon = icon("magic")),
      uiOutput("opttran")
    ),
    mainPanel(
      uiOutput("countsummary"),
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
      uiOutput("countpca"),
      plotlyOutput("pca")
    )
  )
)



# DEG Analysis ----
tab.deg <- tabPanel(
  title = "DEG Analysis ", icon = icon("bar-chart"), 
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("DEG Parameters"),
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
      actionButton("godeg", "Submit", icon = icon("cogs")),
      br(),
      br(),
      br(),
      uiOutput("vishead"),
      uiOutput("vistype"),
      br(),
      uiOutput("deglfcfilt"),
      uiOutput("degpadfilt")
    ),
    mainPanel = mainPanel(
      br(),
      uiOutput("degcomp"),
      plotlyOutput("degplot"),
      br(),
      br(),
      DT::dataTableOutput("mytable"),
      br(),
      br(),
      br()
    )
  )
)



# Help ----
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
