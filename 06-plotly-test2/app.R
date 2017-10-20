#------------------------------------------------------------------------------
# Title:  Shiny Test 06 - Cross Talk
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.19.17
#------------------------------------------------------------------------------

# Package logic ----

## CRAN
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}
packages <- c("shiny", "DT", "tidyr", "plotly")
pack.man(packages)

## Bioconductor
source("https://bioconductor.org/biocLite.R")
if (!require("DESeq2")) biocLite("DESeq2")
if (!require("edgeR")) biocLite("edgeR")



# User interface ----
ui <- fluidPage(
  tabsetPanel(
    tabPanel(
      title = "Welcome", icon = icon("hand-spock-o"),
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
          h4("QC Summary"),
          verbatimTextOutput("contents"),
          br(),
          h4("Count data distributions - box and whisker"),
          plotlyOutput("boxplot"),
          br(),
          h4("Count data distributions - histograms"),
          plotlyOutput("hist"),
          h4("Total reads"),
          plotlyOutput("barplot"),
          br(),
          h4("Principal component analysis"),
          p("Working on it...")
        )
      )
    ),
    tabPanel(
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
          textInput(
            inputId = "degfact",
            label = "Define factor",
            ""
          ),
          textInput(
            inputId = "degfactlev1",
            label = "Factor level 1",
            ""
          ),
          textInput(
            inputId = "degfactlev2",
            label = "Factor level 2",
            ""
          ),
          selectInput(
            inputId = "plottype",
            label = "Choose plot type",
            choices = c(
              "MA plot" = "maplot",
              "Volcano plot" = "volplot"
            )
          ),
          actionButton("godge", "Submit", icon = icon("magic"))          
        ),
        mainPanel = mainPanel(
          br(),
          h4("Comparisons"),
          plotlyOutput("volplot"),
          br(),
          DT::dataTableOutput("mytable"),
          downloadButton("downloadData", "Download Filtered Data")
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
        dds <- dds[ rowSums(counts(dds)) > input$prefilt, ]
        dds <- DESeq(dds)
        
        # reslfc <- lfcShrink(dds = dds, coef = 2, res = results(dds))
        
        return(list(dds))  
      })
    }
  })

  # Pop-up input - transformation method
  output$opttran <- renderUI({
    validate(
      need(input$go != 0, "")
    )
    selectInput(
      inputId = "transform",
      label = HTML(paste("<br/>Choose transformation method for counts")),
      choices = c(
        "Normal log: log2(n + psuedocount)" = "log",
        "Regularized log: rlog(n)" = "rlog",
        "Variance stabilizing transform: vst(n)" = "vst",
        "No transformation" = "raw"
      )       
    )
  })

  # Box plot
  output$boxplot <- renderPlotly({
    validate(
      need(input$go != 0, "")
    )
    if (input$transform == "log") {
      dds <- ddsout()[[1]]
      tmp <- normTransform(dds)
      tmp <- assay(tmp)
      lab <- "log<sub>2</sub>(counts + 1)"
    } else if (input$transform == "rlog") {
      dds <- ddsout()[[1]]
      tmp <- rlog(dds)
      tmp <- assay(tmp)
      lab <- "rlog(counts)"
    } else if (input$transform == "vst") {
      dds <- ddsout()[[1]]
      tmp <- vst(dds)
      tmp <- assay(tmp)
      lab <- "vst(counts)"
    } else if (input$transform == "raw") {
      dds <- ddsout()[[1]]
      tmp <- assay(dds)
      lab <- "Raw counts"
    }
    if (input$go == 0) {
      return()
    } else {
      isolate({
        box <- as.data.frame(tmp)
        box <- tidyr::gather(box)

        plot_ly(
          box,
          type = "box",
          y = ~value,
          color = ~key
        ) %>%
        layout(
          margin = list(b = 90),
          xaxis = list(title = "", tickangle = -45),
          yaxis = list(title = lab)
        )
      })
    }
  })

  output$barplot <- renderPlotly({
    if (input$go == 0) {
      return()
    } else {
      isolate({
        dds <- ddsout()[[1]]
        bar <- as.data.frame(assay(dds))
        bar <- colSums(bar)
        bar <- as.data.frame(t(bar))
        bar <- gather(bar)

        plot_ly(
          bar,
          type = "bar",
          y = ~value,
          x = ~key,
          color = ~key
        ) %>%
        layout(
          margin = list(b = 90),
          xaxis = list(title = "", tickangle = -45),
          yaxis = list(title = "Counts")
        )
      })
    }
  })

  output$hist <- renderPlotly({
    validate(
      need(input$go != 0, "")
    )
    if (input$transform == "log") {
      dds <- ddsout()[[1]]
      tmp <- normTransform(dds)
      tmp <- assay(tmp)
      lab <- "log<sub>2</sub>(counts + 1)"
    } else if (input$transform == "rlog") {
      dds <- ddsout()[[1]]
      tmp <- rlog(dds)
      tmp <- assay(tmp)
      lab <- "rlog(counts)"
    } else if (input$transform == "vst") {
      dds <- ddsout()[[1]]
      tmp <- vst(dds)
      tmp <- assay(tmp)
      lab <- "vst(counts)"
    } else if (input$transform == "raw") {
      dds <- ddsout()[[1]]
      tmp <- assay(dds)
      lab <- "Raw counts"
    }
    if (input$go == 0) {
      return()
    } else {
      isolate({
        hist <- as.data.frame(tmp)
        hist <- tidyr::gather(hist)
        plot_ly(
          data = hist,
          type = "histogram",
          x = ~value,
          color = ~key
        ) %>%
        layout(
          yaxis = list(title = lab)
        )
      })
    }
  })

  # output$optfact1 <- renderUI({
  #   if (input$go == 0) {
  #     return()
  #   } else {
  #     tmp <- ddsout()[[1]]
  #     tmp <- colData(tmp)
  #     validate(
  #       need(input$degfact == names(tmp), "Please enter correct factor")
  #     )
  #     selectInput("degfactlev", "Factor level 1", tmp[, input$degfact])
  #   }
  # })

  # output$optvs <- renderUI({
  #   if (input$go == 0) {
  #     return()
  #   } else {
  #     tmp <- ddsout()[[1]]
  #     tmp <- colData(tmp)
  #     validate(
  #       need(input$degfact == names(tmp), "")
  #     )
  #     h4("Versus")
  #   }
  # })

  # output$optfact2 <- renderUI({
  #   if (input$go == 0) {
  #     return()
  #   } else {
  #     tmp <- ddsout()[[1]]
  #     tmp <- colData(tmp)
  #     validate(
  #       need(input$degfact == names(tmp), "")
  #     )
  #     selectInput("degfactlev", "Factor level 2", tmp[, input$degfact])
  #   }
  # })

  # Volcano plot
  output$volplot <- renderPlotly({
    if (input$godge == 0) {
      return()
    } else {
      isolate({
        if (input$plottype == "maplot") {
        ## Define later
        dds <- ddsout()[[1]]
        fact <- input$degfact
        lev1 <- input$degfactlev1
        lev2 <- input$degfactlev2
        padj <- 0.05
        lfc <- 2
    
    
        ## Automatic data priming
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- subset(datobj, baseMean != 0)
        datobj$isDE <- ifelse(datobj$padj < padj, TRUE, FALSE)
        datobj$isDE[is.na(datobj$isDE)] <- FALSE
        datobj$isLFC <- ifelse(abs(datobj$log2FoldChange) > lfc, TRUE, FALSE)
        datobj$isLFC[is.na(datobj$isLFC)] <- FALSE
    
        ## Tooltips
        tooltips <- paste0(
          "<b>ID:</b> ", rownames(datobj), "<br />",
          "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
          "<b>FDR:</b> ", round(datobj$padj, 3), "<br />",
          "<b>BM:</b> ", round(datobj$baseMean, 3) 
        )    
            
        ## The plot    
        plot_ly(
          type = "scatter",
          mode = "markers",
          x = ~log10(datobj$baseMean),
          y = ~datobj$log2FoldChange,
          color = ~datobj$isDE,
          colors = c("#CCCCCC", "#4286f4"),
          text = tooltips,
          hoverinfo = "text",
          marker = list(size = 3)
        ) %>%
        layout(
          xaxis = list(title = "Base mean"),
          yaxis = list(title = "log<sub>2</sub> fold change")
        )
        # plot_ly(mtcars, x = ~mpg, y = ~wt)
      } else if (input$plottype == "volplot") {
        ## Define later
        dds <- ddsout()[[1]]
        fact <- input$degfact
        lev1 <- input$degfactlev1
        lev2 <- input$degfactlev2
        padj <- 0.05
        lfc <- 2
    
    
        ## Automatic data priming
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- subset(datobj, baseMean != 0)
        datobj$isDE <- ifelse(datobj$padj < padj, TRUE, FALSE)
        datobj$isDE[is.na(datobj$isDE)] <- FALSE
        datobj$isLFC <- ifelse(abs(datobj$log2FoldChange) > lfc, TRUE, FALSE)
        datobj$isLFC[is.na(datobj$isLFC)] <- FALSE
    
        ## Tooltips
        tooltips <- paste0(
          "<b>ID:</b> ", rownames(datobj), "<br />",
          "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
          "<b>FDR:</b> ", round(datobj$padj, 3), "<br />",
          "<b>BM:</b> ", round(datobj$baseMean, 3) 
        )

        plot_ly(
          type = "scatter",
          mode = "markers",
          x = datobj$log2FoldChange,
          y = -log10(datobj$pvalue),
          color = ~datobj$isDE,
          colors = c("#CCCCCC", "#4286f4"),
          text = tooltips,
          hoverinfo = "text",
          marker = list(size = 3)
        ) %>%
        layout(
          xaxis = list(title = "log<sub>2</sub> fold change"),
          yaxis = list(title = "-log<sub>10</sub>(p-value)")
        )
        } 
      }) 
    }
  })
  
  output$mytable <- DT::renderDataTable({
    if (input$godge == 0) {
      return()
    } else {
      isolate({
        dds <- ddsout()[[1]]
        fact <- input$degfact
        lev1 <- input$degfactlev1
        lev2 <- input$degfactlev2
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- round(datobj, 3)
        datobj <- datatable(datobj, filter = "top")
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
}


# Run shiny ----
shinyApp(ui, server)