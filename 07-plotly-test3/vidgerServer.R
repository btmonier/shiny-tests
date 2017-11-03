#------------------------------------------------------------------------------
# Title:  Shiny Test 07 - Modularity - Server Logic
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.26.17
#------------------------------------------------------------------------------

source("server-qc.R")
vidgerServer <- function(input, output) {

  # Load and prime data
  ddsout <- qc(input, output)
  qcopttran(input, output)

  # QC - Header
  output$countsummary <- renderUI({
    if(input$go == 0) {
      return()
    } else {
      h4("QC Summary")
    }
  })

  ## QC - Header - Boxplot
  output$countbox <- renderUI({
    if(input$go == 0) {
      return()
    } else {
      h4("Count data distributions - box and whisker")
    }
  })


  # QC - Boxplot
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

  ## QC - Headers - Histogram
  output$counthist <- renderUI({
    if(input$go == 0) {
      return()
    } else {
      h4("Count data distributions - histogram")
    }
  })

  # QC - Histogram
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

  ## QC - Headers - Barplot
  output$counttotal <- renderUI({
    if(input$go == 0) {
      return()
    } else {
      h4("Total reads")
    }
  })

  # QC - Barplot
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

  ## QC - Headers - PCA
  output$countpca <- renderUI({
    if(input$go == 0) {
      return()
    } else {
      h4("Principal component analysis")
    }
  })


  # ----

  # Conditional levels for DEG
  output$degfact <- renderUI({
    tmp <- input$file2
    tmp <- read.csv(tmp$datapath, header = TRUE, row.names = 1)
    selectInput("fact", "Define factor", colnames(tmp))
  })
  output$deglev1 <- renderUI({
    tmp <- input$file2
    tmp <- read.csv(tmp$datapath, header = TRUE, row.names = 1)
    selectInput("lev1", "Define level 1", tmp[, input$fact])
  })
  output$deglev2 <- renderUI({
    tmp <- input$file2
    tmp <- read.csv(tmp$datapath, header = TRUE, row.names = 1)
    selectInput("lev2", "Define level 2", tmp[, input$fact])
  })

  # Reactive DEG object
  deg <- reactive({
    if (input$go == 0) {
      return()
    } else {
      dds <- ddsout()[[1]]
      # fact <- input$fact
      # lev1 <- input$lev1
      # lev2 <- input$lev2
      dds <- DESeq(dds)
      return(list(dds))
    }
  })

  ## DEG - Headers - Conditional header
  output$degcomp <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
        if (input$plottype == "maplot") {
          h4("Comparisons - MA Plot")
      } else {
        if (input$plottype == "volplot") {
          h4("Comparisons - Volcano Plot")
        }
      }
    } 
  })

  # Crosstalk?
  cross <- SharedData$new(deg()[[1]])

  # DEG plots
  output$degplot <- renderPlotly({
    if (input$godeg == 0) {
      return()
    } else {
      isolate({
        ## Establish objects
        dds <- deg()[[1]]
        fact <- input$fact
        lev1 <- input$lev1
        lev2 <- input$lev2

        ## Automatic data priming
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- subset(datobj, baseMean != 0)
        
        ## Tooltips
        tooltips <- paste0(
          "<b>ID:</b> ", rownames(datobj), "<br />",
          "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
          "<b>FDR:</b> ", round(datobj$padj, 3), "<br />",
          "<b>BM:</b> ", round(datobj$baseMean, 3) 
        )  
    

        if (input$plottype == "maplot") {

          s <- input$mytable_rows_selected  
          if (!length(s)) {  
            p <- cross %>% plot_ly(
              type = "scatter",
              mode = "markers",
              x = ~log10(datobj$baseMean),
              y = ~datobj$log2FoldChange,
              color = I("black"),
              name = "Unfiltered",
              text = tooltips,
              hoverinfo = "text",
              marker = list(size = 3)
            ) %>%
            layout(
              xaxis = list(title = "Base mean"),
              yaxis = list(title = "log<sub>2</sub> fold change")
            ) %>%
            highlight("plotly_selected", color = I("red"), selected = attrs_selected(name = "Filtered"))
          } else if (length(s)) {
            pp <- dds %>% plot_ly() %>%
            add_trace(
              x = ~log10(datobj$baseMean),
              y = ~datobj$log2FoldChange,
              mode = "markers",
              color = I("black")
            ) %>%
            layout(showlegend = TRUE)
            pp <- add_trace(
              pp, 
              data = dds[s, , drop = FALSE],
              x = ~log10(datobj$baseMean),
              y = ~datobj$log2FoldChange,
              mode = "markers",
              color = I("red"),
              name = "Filtered"
            )
          }
        } else if (input$plottype == "volplot") {
          plot_ly(
            type = "scatter",
            mode = "markers",
            x = datobj$log2FoldChange,
            y = -log10(datobj$pvalue),
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
    if (input$godeg == 0) {
      return()
    } else {
      isolate({
        dds <- deg()[[1]]
        fact <- input$fact
        lev1 <- input$lev1
        lev2 <- input$lev2
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- round(datobj, 3)

        datobj2 <- datobj[cross$selection(), ]
        dt <- datatable(datobj, filter = "top")
        if (NROW(datobj2) == 0) {
          dt
        } else {
          formatStyle(
            dt, 
            "rowname", 
            target = "row",
            color = styleEqual(
              datobj2$rowname, rep("white", length(datobj2$rowname))
            ),
            backgroundColor = styleEqual(
              datobj2$rowname, rep("black", length(datobj2$rowname))
            )
          )
        }
        
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