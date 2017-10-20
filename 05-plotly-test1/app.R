# Package logic ----

## CRAN
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("shiny", "DT")
pack.man(packages)

## Bioconductor
source("https://bioconductor.org/biocLite.R")
if (!require("DESeq2")) biocLite("DESeq2")
if (!require("edgeR")) biocLite("edgeR")


# User interface ----
ui <- fluidPage(
  tabsetPanel(
    tabPanel(
      title = "Welcome", icon = icon("bullhorn"),
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          h1("Welcome to the ViDGER2 Pipeline"),
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
    ),
    tabPanel(
      title = "Analyze ", icon = icon("cogs"), 
      fluid = TRUE,
       sidebarLayout(
        sidebarPanel(
          h4("Parameters"),
          selectInput(
            inputId = "dgechoice",
            label = "Choose analysis type",
            choices = c(
              "DESeq2",
              "edgeR",
              "limma"
            ) 
          ),
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
          radioButtons("filt1", "Pre-filter data?", choices = c("Yes", "No")),
          actionButton(inputId = "go", "Submit", icon = icon("magic"))
        ),
        mainPanel(
          br(),
          h4("QC Summary"),
          verbatimTextOutput("contents"),
          br(),
          h4("Comparisons"),
          DT::dataTableOutput("mytable"),
          downloadButton("downloadData", "Download Filtered Data")
        )
      )
    ),
    tabPanel(
      title = "Visualize ", icon = icon("bar-chart"), 
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          h4("Volcano Plot"),
          plotOutput("volPlot", width = "400px", height = "600px"),
          br(),
          h4("MA Plot"),
          plotOutput("maplot", width = "400px", height = "450px")
        )
      )
    ),
    tabPanel(
      title = "Interact ", icon = icon("binoculars"),
      fluid = TRUE,
      sidebarLayout(
        sidebarPanel = NULL,
        mainPanel = mainPanel(
          br(),
          p("Interact with plotly here...")
        )
      )
    ),
    tabPanel(
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

  datasetInput <- reactive({
    ddsout()[[2]]
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
        }  
      })
    }
    summ <- summary(results(dds))
  })  
  
  output$mytable <- DT::renderDataTable({
    if (input$go == 0) {
      return()
    } else {
      isolate({
        reslfc <- as.data.frame(datasetInput())
        reslfc <- round(reslfc, 3)
        reslfc <- datatable(reslfc, filter = "top")
      })  
    }
  })

  output$downloadData <- downloadHandler(
    filename = function() {
      paste("filtered-results", ".csv", sep = "")
    },
    content = function(file) {
      write.csv(
        datasetInput()[input[["mytable_rows_all"]], ], 
        file, row.names = TRUE
      )
    }
  )

  output$volPlot <- renderPlot({
    res <- ddsout()[[1]]
    res <- results(res)
    plot(
      x = res$log2FoldChange, 
      y = -log10(res$pvalue), 
      xlim = c(-2.5, 2.5),
      ylim = c(0, 15),
      xlab = "log2 fold change",
      ylab = "-log10(p value)"
    )
  }) 

  output$maplot <- renderPlot({
    reslfc <- ddsout()[[2]]
    plot(
      x = reslfc$baseMean,
      y = reslfc$log2FoldChange,
      ylim = c(-3, 3),
      xlab = "A",
      ylab = "M"
    )
  }) 
}


# Run shiny ----
shinyApp(ui, server)