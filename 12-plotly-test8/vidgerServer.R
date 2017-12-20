#------------------------------------------------------------------------------
# Title:  Shiny Test 12 - Experimental Designs - Server Logic
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   12.08.17
#------------------------------------------------------------------------------

# Change file upload size to 30 MB
options(shiny.maxRequestSize = 30 * 1024^2)



# Server function
vidgerServer <- function(input, output) {
  
  ## Source exp. method functions
  source("vidger-functions.R")

  ## Example data
  f1 <- as.matrix(read.csv("count-data.csv", header = TRUE, row.names = 1))
  f2 <- read.csv("col-data.csv", header = TRUE, row.names = 1)

  ## Data load option - count data
  output$file1 <- renderUI({
    if (input$examplechoice == "yes") {
      return()
    } else if (input$examplechoice == "no") {
      fileInput(
        inputId = "file1", 
        label = "Submit count data (CSV)",
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      )      
    }
  })

  ## Data load option - metadata
  output$file2 <- renderUI({
    if (input$examplechoice == "yes") {
      return()
    } else if (input$examplechoice == "no") {
      fileInput(
        inputId = "file2", 
        label = "Submit metadata (CSV)",
        accept = c(
          "text/csv",
          "text/comma-separated-values,text/plain",
          ".csv"
        )
      )      
    }
  })

	## QC - reactive - load and add data to DESeqDataSet class
	ddsout <- eventReactive(input$goqc, {
    if (input$examplechoice == "yes") {
      cts <- f1
      coldata <- f2
    } else if (input$examplechoice == "no") {
      cts <- input$file1
      coldata <- input$file2
      cts <- as.matrix(read.csv(cts$datapath, header = TRUE, row.names = 1))
      coldata <- read.csv(coldata$datapath, header = TRUE, row.names = 1)      
    }

    cts <- cts[, rownames(coldata)]
    cts.filt <- cts[rowSums(cts) > input$prefilt, ]
    
    dds <- DESeqDataSetFromMatrix(
      countData = cts.filt,
      colData = coldata,
      design = ~ 1
    )
    return(list(dds, coldata, cts.filt, cts))
	})

  ## QC - reactive - data transformation
  ddstran <- eventReactive(input$goqc, {
    dds <- ddsout()[[1]]
    if (input$transform == "log") {
      tmp <- normTransform(dds)
      lab <- "log<sub>2</sub>(counts + 1)"
    } else if (input$transform == "rlog") {
      tmp <- rlog(dds)
      lab <- "rlog(counts)"
    } else if (input$transform == "vst") {
      tmp <- vst(dds)
      lab <- "vst(counts)"
    } else if (input$transform == "raw") {
      tmp <- dds
      lab <- "Raw counts"
    }
    return(list(tmp, lab))    
  })

  ## QC - header (2) - file summary (count data)
  output$filesummarycts <- renderUI({
    if(input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Count data (first 6 rows)")
    }
  })

  ## QC - header (2) - file summary (col data)
  output$filesummarycoldata <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Sample metadata")
    }
  })

  ## QC - header (2) - ID counts (pre)
  output$headcountpre <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Number of IDs (pre-filtration)")
    }
  })

  ## QC - header (2) - ID counts (post)
  output$headcountpost <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Number of IDs (post-filtration)")
    }
  })

  ### QC - verbatim console - count data head
  output$fileoutputcts <- renderPrint({
    cts <- ddsout()[[3]]
    head(cts)
  })

  ### QC - verbatim console - metadata
  output$fileoutputcoldata <- renderPrint({
    coldata <- ddsout()[[2]]
    coldata
  })

  ### QC - verbatim console - gene no. (pre-filter)
  output$fileoutputcountpre <- renderPrint({
    cts.pre <- ddsout()[[4]]
    nrow(cts.pre)
  })

  ### QC - verbatim console - gene no. (pre-filter)
  output$fileoutputcountpost <- renderPrint({
    cts.post <- ddsout()[[3]]
    nrow(cts.post)
  })

  ## QC - header (2) - boxplot
  output$countbox <- renderUI({
    if(input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Count data distributions - box and whisker")
    }
  })


  ### QC - visualize - boxplot
  output$boxplot <- renderPlotly({
    withProgress(message = "Compiling boxplots...", value = 0, {
      incProgress()
      input$goqc
      tmp <- ddstran()[[1]]
      tmp <- assay(tmp)
      lab <- ddstran()[[2]]
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
    })
  })

  ## QC - header (2) - histogram
  output$counthist <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Count data distributions - histogram")
    }
  })

  ### QC - visualize - histogram
  output$hist <- renderPlotly({
    withProgress(message = "Compiling histograms...", value = 0, {
      incProgress()
      input$goqc
      tmp <- ddstran()[[1]]
      tmp <- assay(tmp)
      lab <- ddstran()[[2]]
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
    })
  })

  ## QC - header (2) - barplot
  output$counttotal <- renderUI({
    if(input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Total reads")
    }
  })

  ### QC - visualize - barplot
  output$barplot <- renderPlotly({
    withProgress(message = "Compiling barplots...", value = 0, {
      incProgress()
      input$goqc
      tmp <- ddstran()[[1]]
      tmp <- assay(tmp)
      lab <- ddstran()[[2]]
      isolate({
        incProgress()
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
    })
  })

  ## QC - header (2) - PCA
  output$headpca <- renderUI({
    if(input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Principal Component Analysis")
    }
  })

  ### QC - select input - choose factor - PCA
  output$pcafact <- renderUI({
    tmp <- ddsout()[[2]]
    selectInput(
      inputId = "pcafact",
      label = "Choose factor",
      choices = colnames(tmp)
    )
  })

  ## QC - visualize - PCA
  output$pca <- renderPlotly({
    tmp <- ddstran()[[1]]
    lab <- ddstran()[[2]]
    validate(
      need(
        expr = class(tmp) == "DESeqTransform",
        message = "Please transform raw counts to visualize PCA."   
      )
    )
    pca.lab <- plotPCA(
      tmp,
      intgroup = input$pcafact
    )
    pca <- plotPCA(
      tmp, 
      intgroup = input$pcafact, 
      returnData = TRUE
    )
    tooltips <- paste0(
      "<b>Sample:</b> ", rownames(pca), "<br />",
      "<b>PC1:</b> ", round(pca$PC1, 3), "<br />",
      "<b>PC2:</b> ", round(pca$PC2, 3)
    )
    plot_ly(
      data = pca,
      type = "scatter",
      mode = "markers",
      x = ~PC1,
      y = ~PC2,
      symbol = ~group,
      marker = list(size = 9),
      text = tooltips,
      hoverinfo = "text"
    ) %>%
    layout(
      xaxis = list(title = pca.lab$labels$x),
      yaxis = list(title = pca.lab$labels$y)
    )   
  })

  ## QC - header (2) - MDS
  output$headmds <- renderUI({
    if(input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Multidimensional Scaling")
    }
  })

  ### QC - select input - choose factor - MDS
  output$mdsfact <- renderUI({
    tmp <- ddsout()[[2]]
    selectInput(
      inputId = "mdsfact",
      label = "Choose factor",
      choices = colnames(tmp)
    )
  })

  ## QC - visualize - MDS
  output$mds <- renderPlotly({
    tmp <- ddstran()[[1]]
    lab <- ddstran()[[2]]
    validate(
      need(
        expr = class(tmp) == "DESeqTransform",
        message = "Please transform raw counts to visualize MDS."   
      )
    )
    mds <- dist(t(assay(tmp)))
    mds <- as.matrix(mds)
    mds <- as.data.frame(colData(tmp)) %>%
      cbind(cmdscale(mds))

    tooltips <- paste0(
      "<b>Sample:</b> ", rownames(mds), "<br />",
      "<b>Coord. 1:</b> ", round(mds[, "1"], 3), "<br />",
      "<b>Coord. 2:</b> ", round(mds[, "2"], 3)
    )
    plot_ly(
      data = mds,
      type = "scatter",
      mode = "markers",
      x = mds[, "1"],
      y = mds[, "2"],
      symbol = mds[, input$mdsfact],
      marker = list(size = 9),
      text = tooltips,
      hoverinfo = "text"
    ) %>%
    layout(
      xaxis = list(title = "MDS coordinate 1"),
      yaxis = list(title = "MDS coordinate 2")
    )   
  })



  ##### ブ レ ー ク  B R E A K  ブ レ ー ク #####



  ## DEG - exp. setup 1 - two group comparisons - factor choice
  output$dgeexp1a <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp1") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp1a",
          label = "Choose factor",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 1 - two group comparisons - choose comparisons
  output$dgeexp1b <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp1") {
        perm <- ddsout()[[2]]
        perm <- as.vector(unique(perm[, input$dgeexp1a]))
        perm <- permutations(n = length(perm), r = 2, v = perm)
        perm <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")
        checkboxGroupInput(
          inputId = "dgeexp1b",
          label = "Choose comparisons you want made",
          choices = as.list(perm)
        )
      }
    }
  })

  ## DEG - exp. setup 2 - mult. group comparisons - factor choice A
  output$dgeexp2a <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
        p(
          em(
            "It appears that your metadata only contains one factor. Please choose 'Two group comparisons' or load data that contains multiple factors."
          ), 
          style = "color:grey"
        )
      } else if (input$dgeexpsetup == "exp2") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp2a",
          label = "Choose factor (A)",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 2 - mult. group comparisons - factor choice B
  output$dgeexp2b <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp2") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp2b",
          label = "Choose factor (B)",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 2 - mult. group comparisons - group comb. choices
  output$dgeexp2c <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp2") {
        perm <- ddsout()[[2]]
        group.c <- factor(
          paste(perm[, input$dgeexp2a], perm[, input$dgeexp2b], sep = "_")
        )
        perm <- levels(group.c)
        perm <- permutations(n = length(perm), r = 2, v = perm)
        perm <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")
        checkboxGroupInput(
          inputId = "dgeexp2c",
          label = "Choose comparisons you want made",
          choices = as.list(perm)
        )
      }
    }
  })

  ## DEG - exp. setup 3 - interaction - factor choice A
  output$dgeexp3a <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
        p(
          em(
            "It appears that your metadata only contains one factor. Please choose 'Two group comparisons' or load data that contains multiple factors."
          ), 
          style = "color:grey"
        )
      } else if (input$dgeexpsetup == "exp3") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp3a",
          label = "Choose factor (A)",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 3 - interaction - factor choice B
  output$dgeexp3b <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp3") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp3b",
          label = "Choose factor (B)",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 3 - interaction - reference level for factor A
  output$dgeexp3c <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp3") {
        tmp <- ddsout()[[2]]
        tmp <- unique(tmp[, input$dgeexp3a])
        selectInput(
          inputId = "dgeexp3c",
          label = "Choose reference level for factor A",
          choices = tmp
        )
      }
    }
  })

  ## DEG - exp. setup 3 - interaction - reference level for factor B
  output$dgeexp3d <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp3") {
        tmp <- ddsout()[[2]]
        tmp <- unique(tmp[, input$dgeexp3b])
        selectInput(
          inputId = "dgeexp3d",
          label = "Choose reference level for factor B",
          choices = tmp
        )
      }
    }
  })

  ## DEG - exp. setup 4 - additive model - blocking factor
  output$dgeexp4a <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
        p(
          em(
            "It appears that your metadata only contains one factor. Please choose 'Two group comparisons' or load data that contains multiple factors."
          ), 
          style = "color:grey"
        )
      } else if (input$dgeexpsetup == "exp4") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp4a",
          label = "Choose your blocking or subjet factor",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 4 - additive model - treatment factor
  output$dgeexp4b <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp4") {
        tmp <- ddsout()[[2]]
        selectInput(
          inputId = "dgeexp4b",
          label = "Choose your treatment factor",
          choices = colnames(tmp)
        )
      }
    }
  })

  ## DEG - exp. setup 4 - additive model - reference level for blocking factor
  output$dgeexp4c <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp4") {
        tmp <- ddsout()[[2]]
        tmp <- unique(tmp[, input$dgeexp4a])
        selectInput(
          inputId = "dgeexp4c",
          label = "Choose reference level for blocking factor",
          choices = tmp
        )
      }
    }
  })

  ## DEG - exp. setup 4 - additive model - reference level for treatment factor
  output$dgeexp4d <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp4") {
        tmp <- ddsout()[[2]]
        tmp <- unique(tmp[, input$dgeexp4b])
        selectInput(
          inputId = "dgeexp4d",
          label = "Choose reference level for treatment factor",
          choices = tmp
        )
      }
    }
  })

  ## DEG - exp. setup - formula - header
  output$dgeexpformhead <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup != "exp1" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp1" & ncol(ddsout()[[2]] < 2)){
        h5(
          strong("Your linear model will look like this:")
        )
      } else {
        h5(
          strong("Your linear model will look like this:")
        )
      }
    }
  })

  ## DEG - exp. setup - formula - formula 1
  output$dgeexpform1 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp1") {
        code(
          paste0(" ~ ", input$dgeexp1a)
        )
      }
    }
  })

  ## DEG - exp. setup - formula - formula 2
  output$dgeexpform2 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp2" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp2") {
        code(
          paste0(" ~ ", input$dgeexp2a, "_", input$dgeexp2b)
        )
      }
    }
  })

  ## DEG - exp. setup 3 - interaction - formula layout
  output$dgeexpform3 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp3" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp3") {
        code(
          paste0(
            " ~ ", input$dgeexp3a,
            " + ", input$dgeexp3b, 
            " + ", input$dgeexp3a,
            ":", input$dgeexp3b
          )
        )
      }
    }
  })

  ## DEG - exp. setup 4 - added effects - formula layout
  output$dgeexpform4 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      if (input$dgeexpsetup == "exp4" & ncol(ddsout()[[2]]) < 2) {
        return()
      } else if (input$dgeexpsetup == "exp4") {
        code(
          paste0(" ~ ", input$dgeexp4a, " + ", input$dgeexp4b)
        )
      }
    }
  })

  ## DEG - edgeR normalization option
  output$dgeexpedgernorm <- renderUI({
    if (input$dgemethod != "edger") {
      return()
    } else {
      selectInput(
        inputId = "dgeexpedgernorm",
        label = "Choose normalization type",
        choices = c(
          "TMM" = "TMM",
          "RLE" = "RLE",
          "upperquartile" = "upperquartile",
          "none" = "none" 
        )
      )
    }
  })

  ## DEG - Choose contrasts
  output$dgemaincontrasts <- renderUI({
    if (input$godge == 0) {
      return()
    } else {
      tmp <- dgeout1()[[2]]
      tmp <- colnames(tmp)
      selectInput(
        inputId = "dgemaincontrasts",
        label = "Choose contrast",
        choices = tmp,
        width = "400px" 
      )
    }
  })

  ## DEG - analysis - reactive
  dgeout1 <- eventReactive(input$godge, {
    cts <- ddsout()[[3]]
    coldata <- ddsout()[[2]]
    if (input$dgemethod == "limma") {
      if (input$dgeexpsetup == "exp1") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          de.genes <- limma.exp1(
            fact = input$dgeexp1a,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp1b
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp2") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          de.genes <- limma.exp2(
            fact1 = input$dgeexp2a,
            fact2 = input$dgeexp2b,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp2c
          )
          fit.names <- de.genes[[2]]       
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp3") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)     
          de.genes <- limma.exp3(
            fact1 = input$dgeexp3a,
            fact2 = input$dgeexp3b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp3c,
            fact2.rlvl = input$dgeexp3d
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp4") {
        withProgress(message = "Running limma-voom...", value = 0, {
          incProgress(1/2)
          de.genes <- limma.exp4(
            fact1 = input$dgeexp4a,
            fact2 = input$dgeexp4b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp4c,
            fact2.rlvl = input$dgeexp4d
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      }
    } else if (input$dgemethod == "edger") {
      if (input$dgeexpsetup == "exp1") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp1(
            fact = input$dgeexp1a,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp1b,
            norm = input$dgeexpedgernorm
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp2") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp2(
            fact1 = input$dgeexp2a,
            fact2 = input$dgeexp2b,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp2c,
            norm = input$dgeexpedgernorm
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp3") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp3(
            fact1 = input$dgeexp3a,
            fact2 = input$dgeexp3b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp3c,
            fact2.rlvl = input$dgeexp3d,
            norm = input$dgeexpedgernorm
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp4") {
        withProgress(message = "Running edgeR...", value = 0, {
          incProgress(1/2)
          de.genes <- edger.exp4(
            fact1 = input$dgeexp4a,
            fact2 = input$dgeexp4b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp4c,
            fact2.rlvl = input$dgeexp4d,
            norm = input$dgeexpedgernorm
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      }
    } else if (input$dgemethod == "deseq") {
      if (input$dgeexpsetup == "exp1") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          de.genes <- deseq.exp1(
            fact = input$dgeexp1a,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp1b 
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp2") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          de.genes <- deseq.exp2(
            fact1 = input$dgeexp2a,
            fact2 = input$dgeexp2b,
            cts = cts,
            coldata = coldata,
            perm.h = input$dgeexp2c
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp3") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          de.genes <- deseq.exp3(
            fact1 = input$dgeexp3a,
            fact2 = input$dgeexp3b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp3c,
            fact2.rlvl = input$dgeexp3d        
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      } else if (input$dgeexpsetup == "exp4") {
        withProgress(message = "Running DESeq2...", value = 0, {
          incProgress(1/2)
          de.genes <- deseq.exp4(
            fact1 = input$dgeexp4a,
            fact2 = input$dgeexp4b,
            cts = cts,
            coldata = coldata,
            fact1.rlvl = input$dgeexp4c,
            fact2.rlvl = input$dgeexp4d
          )
          fit.names <- de.genes[[2]]
          fit.cont <- de.genes[[1]]
          incProgress(2/2)
        })
      }
    } 
    return(list(fit.cont, fit.names))
  })

  ## DGE - header (2) - DGE Overview
  output$headdgeoverview <- renderUI({
    if(input$godge == 0) {
      p(
        br(),
        em(
          "Choose an experimental setup with parameters and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("DGE Overview")
    }
  })

  ## DGE - header (2) - DGE Overview
  output$headdgeplots <- renderUI({
    if(input$godge == 0) {
      p(
        br(),
        em(
          "Choose an experimental setup with parameters and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Interactive Plots")
    }
  }) 

  ## DGE - contrast table
  dgeout2 <- reactive({
    if(input$godge == 0) {
      return()
    } else {
      isolate({
        expset <- input$dgeexpsetup
      })
      contTable <- getContTable(
        de.genes = dgeout1()[[1]],
        coef = input$dgemaincontrasts,
        cts = ddsout()[[3]],
        expset = expset,
        design = dgeout1()[[2]],
        fact = input$dgeexp1a
      )
      return(list(contTable))
    }
  })  

  ## DEG - select input - Visualization type
  output$vistype <- renderUI({
    if(input$godge == 0) {
      return()
    } else {
      radioButtons(
        inputId = "plottype",
        label = "Choose plot type",
        choices = c(
          "MA plot" = "maplot",
          "Volcano plot" = "volplot"
        ),
        selected = "maplot",
        inline = TRUE
      )
    }
  })
  
  # DEG - Filter Data 
  dgeout3 <- reactive({
    if (input$godge == 0) {
      return()
    } else {
      tmp <- dgeout2()[[1]]
      tmp$padj <- round(tmp$padj, 3)
      tmp[is.na(tmp)] <- 1
      padj <- input$dgepadjcutoff
      lfc <- input$dgefcmin
      tmp <- tmp[abs(tmp$log2FoldChange) >= lfc, ]
      tmp <- tmp[tmp$padj <= padj, ]
      return(tmp)
    }
  })

  ## DGE - DEBUG
  output$debugdge <- renderPrint({
    dim(dgeout3())
  })

  # DEG - Shared data
  share <- reactive({
    SharedData$new(dgeout3())
  })

  # DEG - Conditional plot
  output$dgeplot <- renderPlotly({
    if (input$godge == 0) {
      return()
    } else {
      s <- input$mytable_rows_selected
      datobj <- dgeout3()
      tooltips <- paste0(
        "<b>ID:</b> ", datobj$id, "<br />",
        "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
        "<b>PADJ:</b> ", round(datobj$padj, 3), "<br />",
        "<b>BM:</b> ", round(datobj$baseMean, 3) 
      )
      if (input$godge == 0) {
        return()
      } else {
        if (input$plottype == "maplot") {
          if (!length(s)) {
            p <- share() %>%
              plot_ly(
                x = ~log10(baseMean), 
                y = ~log2FoldChange,
                type = "scatter", 
                mode = "markers", 
                color = I("darkgray"),
                text = tooltips,
                hoverinfo = "text",
                name = "Unfiltered") %>%
              layout(
                showlegend = TRUE,
                xaxis = list(title = "log<sub>10</sub>(baseMean)"),
                yaxis = list(title = "log<sub>2</sub>(fold change)")
              ) %>% 
              highlight(
                "plotly_selected", 
                color = I("royalblue1"), 
                selected = attrs_selected(
                  name = "Filtered"
                )
              )
          } else if (length(s)) {
            pp <- dgeout3() %>%
              plot_ly() %>% 
              add_trace(
                x = ~log10(baseMean), 
                y = ~log2FoldChange,
                type = "scatter", 
                mode = "markers", 
                color = I("darkgray"),
                text = tooltips,
                hoverinfo = "text", 
                name = "Unfiltered") %>%
              layout(
                showlegend = TRUE,
                xaxis = list(title = "log<sub>10</sub>(baseMean)"),
                yaxis = list(title = "log<sub>2</sub>(fold change)")
              )
    
            # selected data
            pp <- add_trace(
              pp, 
              data = dgeout3()[s, , drop = FALSE], 
              x = ~log10(baseMean), 
              y = ~log2FoldChange,
              type = "scatter", 
              mode = "markers", 
              color = I("royalblue1"), 
              # text = tooltips,
              size = I(8), 
              name = "Filtered"
            )
          }
        } else if (input$plottype == "volplot") {
          if (!length(s)) {
            p <- share() %>%
              plot_ly(
                x = ~log2FoldChange, 
                y = ~-log10(pvalue),
                type = "scatter", 
                mode = "markers", 
                color = I("darkgray"),
                text = tooltips,
                hoverinfo = "text", 
                name = "Unfiltered") %>%
              layout(
                showlegend = TRUE,
                xaxis = list(title = "log<sub>2</sub>(fold change)"),
                yaxis = list(title = "-log<sub>10</sub>(p-value)")
              ) %>%
              highlight(
                "plotly_selected", 
                color = I("royalblue1"), 
                selected = attrs_selected(
                  name = "Filtered"
                )
              )
          } else if (length(s)) {
            pp <- dgeout3() %>%
              plot_ly() %>% 
              add_trace(
                x = ~log2FoldChange, 
                y = ~-log10(pvalue),
                type = "scatter", 
                mode = "markers", 
                color = I("darkgray"),
                text = tooltips,
                hoverinfo = "text", 
                name = "Unfiltered") %>%
              layout(
                showlegend = TRUE,
                xaxis = list(title = "log<sub>2</sub>(fold change)"),
                yaxis = list(title = "-log<sub>10</sub>(p-value)")
              )
    
            # selected data
            pp <- add_trace(
              pp, 
              data = dgeout3()[s, , drop = FALSE], 
              x = ~log2FoldChange, 
              y = ~-log10(pvalue), 
              mode = "markers", 
              color = I("royalblue1"),
              # text = tooltips, 
              size = I(8), 
              name = "Filtered"
            )
          }      
        }
      }
    }
  })


  # DEG - Conditional table
  output$mytable <- DT::renderDataTable({
    if (input$godge == 0) {
      return()
    } else {
      tmp <- dgeout3() %>% mutate_if(is.numeric, round, digits = 3)
      m2 <- tmp[share()$selection(),]
      dt <- DT::datatable(tmp)
      if (NROW(m2) == 0) {
        dt
      } else {
        DT::formatStyle(
          dt, 
          "id", 
          target = "row", 
          color = DT::styleEqual(m2$id, rep("white", length(m2$id))),
          backgroundColor = DT::styleEqual(
            m2$id, 
            rep("darkgray", length(m2$id))
          )
        )
      }
    }
  }) 

  # DEG - Show download button (FILTERED)
  output$downloadfilt <- renderUI({
    validate(
      need(
        expr = !is.null(input$godge),
        message = ""   
      )
    )
    if(input$godge == 0) {
      return()
    } else {
      downloadButton("downfiltData", "Download Filtered Data")
    }     
  })

  # Download FILTERED data
  output$downfiltData <- downloadHandler(
    filename = function() {
      paste0(
        input$dgemaincontrasts,
        "filtered_results_", 
        nrow(dgeout3()[input[["mytable_rows_all"]], ]),
        ".csv"
      )
    },
    content = function(file) {
      write.csv(
        dgeout3()[input[["mytable_rows_all"]], ], 
        file, row.names = FALSE
      )
    }
  )   

  # DEG - Show download button (ALL)
  output$downloadall <- renderUI({
    validate(
      need(
        expr = !is.null(input$godge),
        message = ""   
      )
    )
    if(input$godge == 0) {
      return()
    } else {
      downloadButton("downallData", "Download All Data")
    }     
  })

  # Download ALL data
  output$downallData <- downloadHandler(
    filename = function() {
      paste0(input$dgemaincontrasts, "_all_results.csv")
    },
    content = function(file) {
      write.csv(
        dgeout2()[[1]], 
        file, row.names = FALSE
      )
    }
  )   



  ### ブ レ ー ク  B R E A K  ブ レ ー ク ###



  ## HEAT - header (2) - Heatmap analysis
  output$headheat <- renderUI({
    if(input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Interactive Heatmap")
    }
  })


  ## HEAT - input - Choose number of variable IDs
  output$heatnumber <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      textInput(
        inputId = "heatnumber",
        label = "ID cutoff",
        value = 20
      )
    }
  })

  ## HEAT - data - choose transformation
  heattran <- reactive({
    if (input$goqc == 0) {
      return()
    } else {
      dds <- ddsout()[[1]]
      if (input$transform == "log") {
        tmp <- normTransform(dds)
        tmp <- assay(tmp)
        lab <- "log<sub>2</sub>(counts + 1)"
      } else if (input$ransform == "rlog") {
        tmp <- rlog(dds)
        tmp <- assay(tmp)
        lab <- "rlog(counts)"
      } else if (input$transform == "vst") {
        tmp <- vst(dds)
        tmp <- assay(tmp)
        lab <- "vst(counts)"
      } else if (input$transform == "raw") {
        tmp <- dds
        tmp <- counts(dds)
        lab <- "Raw counts"
      }
    }
    return(list(tmp, lab))    
  })

  ## HEAT - data - get variable IDs
  heattran2 <- reactive({
    if (input$goqc == 0) {
      return()
    } else {
      dds.counts <- ddsout()[[3]]
      heat.counts <- heattran()[[1]]
      num <- input$heatnumber
      topID <- order(rowMeans(dds.counts), decreasing = TRUE)
      heat.mat <- heat.counts[topID, ]
      # heat.mat <- heat.mat - rowMeans(heat.mat)
      heat.mat <- heat.mat[1:num, ,drop = FALSE]
    }
    return(list(heat.mat))
  })




  ## HEAT - visualization - Plotly heatmap
  output$heatplot1 <- renderPlotly({
    if (input$goqc == 0) {
      return()
    } else {
      heat <- heattran2()[[1]]
      tooltips <- paste0(
        # "<b>Sample:</b> ", colnames(heat), "<br />",
        "<b>ID:</b> ", rownames(heat), "<br />",
        "<b>Count:</b> ", round(heat, 3)
      )
      tooltips <- matrix(tooltips, ncol = ncol(heat), byrow = FALSE)
      plot_ly(
        x = colnames(heat),
        y = rownames(heat),
        z = heat,
        type = "heatmap",
        text = tooltips,
        hoverinfo = "text",
        source = "heatplot",
      ) %>%
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        margin = list(l = 100)
      )
    }
  })

  ### HEAT - select input - choose factor - PCA
  output$heatfactor <- renderUI({
    tmp <- ddsout()[[2]]
    selectInput(
      inputId = "heatfactor",
      label = "Choose factor",
      choices = colnames(tmp)
    )
  })

  ## HEAT - visualization - count plot
  output$heatplot2 <- renderPlotly({
    if (input$goqc == 0) {
      return()
    } else {
      s <- event_data("plotly_click", source = "heatplot")
      # source("misc-functions.R")
      rc.data <- counts(ddsout()[[1]])
      test <- getGenes(
        rc.data = rc.data, 
        id = s[["y"]],
        coldata = ddsout()[[2]]
      )
      tooltips <- paste0(
        "<b>Sample:</b> ", test$sample, "<br />",
        "<b>Counts:</b> ", round(test$counts, 3)
      )

      plot_ly(
        data = test,
        type = "scatter",
        mode = "markers",
        x = test[, input$heatfactor],
        y = test[, "counts"],
        color = test[, input$heatfactor],
        text = tooltips,
        marker = list(size = 7),
        hoverinfo = "text" 
      ) %>%
      layout(
        title = paste(s[["y"]], "Counts"),
        xaxis = list(title = paste(input$fact)),
        yaxis = list(title = "Normalized counts")
      )
    }
  })



  ### B R E A K ###



  ## BIC - header (2) - Bicluster analysis
  output$headbic <- renderUI({
    if (input$goqc == 0) {
      p(
        br(),
        em(
          "Load data and click the 'submit' button to see the results."
        ), 
        style = "color:grey"
      )
    } else {
      h4("Bicluster Analysis")
    }
  })

  ## BIC - header (3) - Bicluster analysis parameters
  output$headbicparameters <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      h5(strong("Parameters"))
    }
  })

  ## BIC - input - Choose number of variable genes
  output$bicvarnumber <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      textInput(
        inputId = "bicvarnumber",
        label = "Variable cutoff value",
        value = 500
      )
    }
  })

  ## BIC - input - Choose bicluster algorithm
  output$bicalg <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      selectInput(
        inputId = "bicalg",
        label = "Choose bicluster algorithm",
        choices = c(
          "QUBIC" = "qubic",
          "Bimax" = "bimax",
          "CC" = "cc",
          "Plaid" = "plaid",
          "Spectral" = "spectral",
          "Xmotifs" = "xmotifs"
        )
      )
    }
  })

  ## BIC - actionbutton - submit biclustering
  output$gobic <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      actionButton("gobic", "Launch Analysis", icon = icon("space-shuttle"))
    }
  })

  ## BIC - reactive - get variable counts
  bicout <- eventReactive(input$gobic, {
    num <- input$bicvarnumber
    cts <- ddsout()[[1]]
    cts <- assay(cts)
    tran <- ddstran()[[1]]
    tran <- assay(tran)
    topID <- order(rowVars(cts), decreasing = TRUE)
    cts.var <- tran[topID, ]
    cts.var <- cts.var[1:num, ]
    if (input$bicalg == "qubic") {
      withProgress(message = "Running BCQU...", value = 0, {
        incProgress(1/2)
        res <- biclust::biclust(x = cts.var, method = BCQU() )
        incProgress(2/2)
      })
    } else if (input$bicalg == "bimax") {
      withProgress(message = "Running BCBimax...", value = 0, {
        incProgress(1/2)
        res <- biclust::biclust(x = cts.var, method = BCBimax() )
        incProgress(2/2)
      })
    } else if (input$bicalg == "cc") {
      withProgress(message = "Running BCCC...", value = 0, {
        incProgress(1/2)
        res <- biclust::biclust(x = cts.var, method = BCCC() )
        incProgress(2/2)
      })
    } else if (input$bicalg == "plaid") {
      withProgress(message = "Running BCPlaid...", value = 0, {
        incProgress(1/2)
        res <- biclust::biclust(x = cts.var, method = BCPlaid() )
        incProgress(2/2)
      })
    } else if (input$bicalg == "spectral") {
      withProgress(message = "Running BCSpectal...", value = 0, {
        incProgress(1/2)
        res <- biclust::biclust(x = cts.var, method = BCSpectral() )
        incProgress(2/2)
      })
    } else if (input$bicalg == "xmotifs") {
      withProgress(message = "Running BCXmotifs...", value = 0, {
        incProgress(1/2)
        res <- biclust::biclust(x = cts.var, method = BCXmotifs() )
        incProgress(2/2)
      })
    }            
    return(list(res, cts.var))
  })

  ## BIC - input - choose cluster
  output$bicclustnumber <- renderUI({
    res <- bicout()[[1]]
    if (res@Number < 1) {
      return()
    } else {
      n <- res@Number
      selectInput(
        inputId = "bicclustnumber",
        label = "Choose cluster",
        choices = c(1:n)
      )
    }
  })

  ### BIC - head (3) - Clust overview
  output$headbicsummary <- renderUI({
    validate(
      need(
        input$bicclustnumber != "",
        ""
      )
    )
    h5(strong("Clustering Overview"))
  })


  ## BIC - text - cluster overview
  output$bicsummary <- renderUI({
    res <- bicout()[[1]]
    if (res@Number < 1) {
      p("No clusters found using this algorithm!")
    } else {
      rn <- res@RowxNumber
      cn <- res@NumberxCol
      clust <- input$bicclustnumber
      clust <- as.numeric(clust)
      ids <- length(rn[, clust][rn[, clust] == TRUE])
      samp <- length(cn[clust, ][cn[clust, ] == TRUE])

      out <- paste0(
        "This algorithm found ", ids, " IDs amongst ", samp, 
        " samples in cluster ", clust, "." 
      )
      p(paste(out))
    }
  })

  ### BIC - head (3) - Clust overview
  output$headbicsummary <- renderUI({
    validate(
      need(
        input$bicclustnumber != "",
        ""
      )
    )
    h5(strong("Clustering Overview"))
  })

  ## BIC - header (3) - Bicluster analysis parameters
  output$headbicheatsummary <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      res <- bicout()[[1]]
      if (res@Number < 1) {
        return()
      } else {
        clust <- input$bicclustnumber
        h5(strong(paste0("Heatmap analysis of cluster ", clust)))
      }
    }
  })

  ## BIC - visualize - Bicluster heatmap
  output$bicheatplot <- renderPlot({
    res <- bicout()[[1]]
    if (res@Number < 1) {
      return()
    } else {
      n <- input$bicclustnumber
      n <- as.numeric(n)
      res <- bicout()[[1]]
      cts.var <- bicout()[[2]]
      par(mar = c(10, 6, 3, 5) + 0.1)
      quheatmap(
        x = cts.var,
        bicResult = res,
        number = n, 
        showlabel = TRUE
      )
    }
    
  })

  # output$debugdge <- renderPrint({
  #   if (input$goqc == 0) {
  #     return()
  #   } else {
  #     tmp <- bicout()[[2]]
  #     summary(tmp)
  #   }
  # })
}