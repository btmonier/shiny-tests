#------------------------------------------------------------------------------
# Title:  Shiny Test 11 - Crosstalk - Tabs
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.24.17
#------------------------------------------------------------------------------

# Welcome page ----
tab.welcome <- tabPanel(
  title = "Welcome", icon = icon("hand-spock-o"),
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      HTML("<h1>Welcome to the ViDGER<sup>2</sup> DGE Pipeline</h1>"),
      h3("About"),
      HTML("<p>ViDGER<sup>2</sup> (<b>Vi</b>sualization of <b>D</b>ifferential <b>G</b>ene <b>E</b>xpression using <b>R</b>), is a web-based tool for the analysis of RNA-seq count data. ViDGER<sup>2</sup> is derived from the R package, ViDGER, which can produce information-rich visualizations for the interpreation of differential gene expression (DGE) results from three widely-used tools, <i>Cuffdiff</i>, <i>DESeq2</i>, and <i>edgeR</i>. ViDGER<sup>2</sup> is a <b>user-friendly</b> and <b>interactive</b> Shiny app for gene expression analysis. This app takes advantage of several popular DGE tools (<i>DESeq2</i>, <i>edgeR</i>, and <i>limma</i>) available through Bioconductor in conjunction with the Plotly and DataTable API libraries for R. </p>"),
      h3("Submit and QC"),
      p("To use this app, all you need to upload is two CSV files: a file of raw read counts and another file for sample treatment data. For further information of how these files should look like, take a look at the demonstration data or go to the help section for a detailed walkthrough."),
      h3("DGE Analysis"),
      p("To determine DGE, you currently have the option of using 3 Bioconductor packages mentioned in the prior sections. DGE can be visualized using common techniques, including volcano plots, MA plots, and heatmaps. Additionally, DGE data sets can be downloaded as CSV files for further analysis in this section."),
      h3("Help"),
      p("Need a walkthrough? Look through this tab for detailed information about the pipeline."),
      h3("FAQ"),
      p("Need a quick reminder of some arbitrary concepts? Want to know what you can and can't do with this app? Head on over to the FAQ tab for some frequently asked question.")
    )
  )
)



# Submit and QC page ----
tab.submit <- tabPanel(
  title = "Submit and QC ", icon = icon("cogs"), 
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("Submission Parameters"),
      radioButtons(
        inputId = "examplechoice",
        label = "How do you want to start?",
        choices = c(
          "Start with some example data." = "yes",
          "Load my own data." = "no"  
        )
      ),
      uiOutput("file1"),
      uiOutput("file2"),
      # fileInput(
      #   inputId = "file1", 
      #   label = "Submit count data (CSV)",
      #   accept = c(
      #     "text/csv",
      #     "text/comma-separated-values,text/plain",
      #     ".csv"
      #   )
      # ),
      # fileInput(
      #   inputId = "file2", 
      #   label = "Submit sample information (CSV)",
      #   accept = c(
      #     "text/csv",
      #     "text/comma-separated-values,text/plain",
      #     ".csv"
      #   )
      # ),
      textInput(
        inputId = "prefilt",
        labe = "Filter cutoff",
        value = ""
      ),
      br(),
      selectInput(
        inputId = "transform",
        label = "Choose transformation method for counts",
        choices = c(
          "Normal log: log2(n + pseudocount)" = "log",
          "Regularized log: rlog(n)" = "rlog",
          "Variance stabilizing transform: vst(n)" = "vst",
          "No transformation" = "raw"
        )       
      ),
      actionButton(inputId = "goqc", "Submit", icon = icon("magic")),
      uiOutput("opttran")
    ),
    mainPanel(
      uiOutput("filesummary"),
      uiOutput("filesummarycts"),
      verbatimTextOutput("fileoutputcts"),
      uiOutput("filesummarycoldata"),
      verbatimTextOutput("fileoutputcoldata"),
      br(),
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
      uiOutput("pcafact"),
      plotlyOutput("pca"),
      br(),
      br(),
      br()
    )
  )
)



# DEG Analysis ----
tab.deg <- tabPanel(
  title = "DGE Analysis 1 ", icon = icon("bar-chart"), 
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("DGE Parameters"),
      selectInput(
        inputId = "dgetype",
        label = "Choose analysis type",
        choices = c(
          "DESeq2" = "deseq",
          "edgeR" = "edger",
          "limma-voom" = "limma"
        )
      ),
      uiOutput("degfact"),
      uiOutput("deglev1"),
      uiOutput("deglev2"),
      uiOutput("deseqform"),
      actionButton("godeg", "Submit", icon = icon("magic")),
      br(),
      br(),
      uiOutput("degvishead"),
      uiOutput("vistype"),
      uiOutput("deglfcfilt"),
      uiOutput("degpadfilt"),
      uiOutput("govisfilt")
    ),
    mainPanel = mainPanel(
      uiOutput("degcomp"),
      uiOutput("degcheck"),
      verbatimTextOutput("debugdeg"),
      plotlyOutput("degplot"),
      br(),
      br(),
      DT::dataTableOutput("mytable"),
      div(style = "display:inline-block", uiOutput("downloadfilt")),
      div(style = "display:inline-block", uiOutput("downloadall")),
      br(),
      br(),
      br()
    )
  )
)



# Heatmap
tab.heat <- tabPanel(
  title = "DGE Analysis 2 ",
  icon = icon("fire"),
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("Heatmap Parameters"),
      textInput(
        inputId = "heatrows",
        label = "Choose cutoff for most variable IDs",
        value = 10
      ),
      selectInput(
        inputId = "heattransform",
        label = "Choose transformation method for counts",
        choices = c(
          "Normal log: log2(n + pseudocount)" = "log",
          "Regularized log: rlog(n)" = "rlog",
          "Variance stabilizing transform: vst(n)" = "vst",
          "No transformation" = "raw"
        )       
      ),
      actionButton(inputId = "goheat", "Submit", icon = icon("magic"))
    ),
    mainPanel = mainPanel(
      uiOutput("heattitle"),
      plotlyOutput("heatplot"),
      br(),
      br(),
      plotlyOutput("heatcount"),
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
      p("Help placeholder")
    )
  )
)



# FAQ ----
tab.faq <- tabPanel(
  title = "FAQ ",
  icon = icon("question-circle"),
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      br(),
      p("FAQ placeholder")
    )
  )
)