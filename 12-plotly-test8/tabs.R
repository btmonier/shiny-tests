#------------------------------------------------------------------------------
# Title:  Shiny Test 12 - Experimental Designs - Tabs
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   12.08.17
#------------------------------------------------------------------------------

# Welcome page ----
tab.welcome <- tabPanel(
  title = "Welcome", icon = icon("hand-spock-o"),
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel = NULL,
    mainPanel = mainPanel(
      HTML("<h1>Welcome to the ViDGER Pipeline</h1>"),
      h3("About"),
      HTML("<p>ViDGER (<b>Vi</b>sualization of <b>D</b>ifferential <b>G</b>ene <b>E</b>xpression using <b>R</b>), is a web-based tool for the analysis of RNA-seq count data. ViDGER is derived from the R package, ViDGER, which can produce information-rich visualizations for the interpreation of differential gene expression (DGE) results from three widely-used tools, <i>Cuffdiff</i>, <i>DESeq2</i>, and <i>edgeR</i>. ViDGER is a <b>user-friendly</b> and <b>interactive</b> Shiny app for gene expression analysis. This app takes advantage of several popular DGE tools (<i>DESeq2</i>, <i>edgeR</i>, and <i>limma</i>) available through Bioconductor in conjunction with the Plotly and DataTable API libraries for R. </p>"),
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
      h4("1. Submission Parameters"),
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
      br(),
      h4("2. Data Processing"),
      textInput(
        inputId = "prefilt",
        label = withMathJax("Filter cutoff (count data row sums \\( < n\\))"),
        value = 10
      ),
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
      br(),
      h4("3. Launch Overview"),
      actionButton(inputId = "goqc", "Submit", icon = icon("space-shuttle")),
      br(),
      p(
        "After you click 'submit', you may proceed to the 'DEG Analysis' tab for differential gene expression analyses."
      )
    ),
    mainPanel(
      tabsetPanel(
        tabPanel(
          title = "File Summary",
          uiOutput("filesummarycts"),
          verbatimTextOutput("fileoutputcts"),
          uiOutput("filesummarycoldata"),
          verbatimTextOutput("fileoutputcoldata"),
          br(),
          uiOutput("headcountpre"),
          verbatimTextOutput("fileoutputcountpre"),
          uiOutput("headcountpost"),
          verbatimTextOutput("fileoutputcountpost"),
          br(),
          br()
        ),
        tabPanel(
          title = "Count Summary",
          uiOutput("countbox"),
          plotlyOutput("boxplot"),
          br(),
          uiOutput("counthist"),
          plotlyOutput("hist"),
          br(),
          uiOutput("counttotal"),
          plotlyOutput("barplot"),
          br(),
          br()
        ),
        tabPanel(
          title = "Correlation",
          uiOutput("headcor"),
          # verbatimTextOutput("debugdge"),
          plotlyOutput("corplot1"),
          br(),
          plotlyOutput("corplot2")
        ),
        tabPanel(
          title = "PCA",
          uiOutput("headpca"),
          uiOutput("pcafact"),
          plotlyOutput("pca"),
          br(),
          br()
        ),
        tabPanel(
          title = "MDS",
          uiOutput("headmds"),
          uiOutput("mdsfact"),
          plotlyOutput("mds"),
          br(),
          br()
        ),
        tabPanel(
          title = "Heatmap",
          uiOutput("headheat"),
          uiOutput("heatnumber"),
          plotlyOutput("heatplot1"),
          br(),
          uiOutput("heatfactor"),
          plotlyOutput("heatplot2")
        ),
        tabPanel(
          title = "Biclustering",
          uiOutput("headbic"),
          uiOutput("headbicparameters"),
          uiOutput("bicvarnumber"),
          uiOutput("bicalg"),
          uiOutput("gobic"),
          br(),
          br(),
          uiOutput("headbicsummary"),
          uiOutput("bicclustnumber"),
          uiOutput("bicsummary"),
          br(),
          uiOutput("headbicheatsummary"),
          plotOutput("bicheatplot", height = 800, width = 600),
          div(style = "display:inline-block", uiOutput("downloadbicfilt")),
          div(style = "display:inline-block", uiOutput("downloadbicplotpdf")),
          br(),
          br(),
          br()
        )  
      )
    )
  )
)



# DEG Analysis ----
tab.deg <- tabPanel(
  title = "DEG Analysis ", icon = icon("bar-chart"), 
  fluid = TRUE,
  sidebarLayout(
    sidebarPanel(
      h4("1. Experimental Setup"),
      selectInput(
        inputId = "dgeexpsetup",
        label = "Choose an experimental design",
        choices = c(
          "Two group comparisons" = "exp1",
          "Multiple factor comparisons (factorial)" = "exp2",
          "Classical interaction design" = "exp3",
          "Additive models - paired or blocking" = "exp4"
        ),
        selected = "",
        multiple = FALSE
      ),
      uiOutput("dgeexp1a"),
      uiOutput("dgeexp1b"),
      uiOutput("dgeexp2a"),
      uiOutput("dgeexp2b"),
      uiOutput("dgeexp2c"),
      uiOutput("dgeexp3a"),
      uiOutput("dgeexp3b"),
      uiOutput("dgeexp3c"),
      uiOutput("dgeexp3d"),
      uiOutput("dgeexp4a"),
      uiOutput("dgeexp4b"),
      uiOutput("dgeexp4c"),
      uiOutput("dgeexp4d"),
      uiOutput("dgeexpformhead"),
      uiOutput("dgeexpform1"),
      uiOutput("dgeexpform2"),
      uiOutput("dgeexpform3"),
      uiOutput("dgeexpform4"),
      br(),
      h4("2. DGE Parameters"),
      selectInput(
        inputId = "dgemethod",
        label = "Choose method",
        choices = c(
          "DESeq2" = "deseq",
          "edgeR" = "edger",
          "limma-voom" = "limma"
        )
      ),
      splitLayout(
        textInput(
          inputId = "dgepadjcutoff",
          label = withMathJax("Adj. \\(p\\)-value cutoff"),
          value = 0.05
        ),
        textInput(
          inputId = "dgefcmin",
          label = "Min. fold change",
          value = 1
        )
      ),
      uiOutput("dgeexpedgernorm"),
      br(),
      h4("3. Launch Analysis"),
      actionButton("godge", "Submit", icon = icon("space-shuttle"))
    ),
    mainPanel = mainPanel(
      tabsetPanel(
        tabPanel(
          title = "Overview",
          uiOutput("headdgeoverview")
        ),
        tabPanel(
          title = "Plots",
          uiOutput("headdgeplots"),
          uiOutput("dgemaincontrasts"),
          uiOutput("vistype"),
          plotlyOutput("dgeplot"),
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