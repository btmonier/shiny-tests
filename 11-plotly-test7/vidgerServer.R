#------------------------------------------------------------------------------
# Title:  Shiny Test 11 - Crosstalk - Server Logic
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.24.17
#------------------------------------------------------------------------------

vidgerServer <- function(input, output) {

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
 


	# QC - Reactive
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
    cts <- cts[rowSums(cts) > input$prefilt, ]
    
    dds <- DESeqDataSetFromMatrix(
      countData = cts,
      colData = coldata,
      design = ~ 1
    )
    return(list(dds, coldata, cts))
	})

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

  # QC - Header - File Summary
  output$filesummary <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h2("File Summary")
    }
  })

  # QC - Header - File Summary (count data)
  output$filesummarycts <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Count data (first 6 rows)")
    }
  })

  # QC - Header - File Summary (col data)
  output$filesummarycoldata <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Sample metadata")
    }
  })


  # QC - verbatim - file summary
  output$fileoutputcts <- renderPrint({
    cts <- ddsout()[[3]]
    head(cts)
  })

  # QC - verbatim - coldata summary
  output$fileoutputcoldata <- renderPrint({
    coldata <- ddsout()[[2]]
    coldata
  })


  # QC - Header - QC Summary
  output$countsummary <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h2("QC Summary")
    }
  })

  ## QC - Header - Boxplot
  output$countbox <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Count data distributions - box and whisker")
    }
  })


  # QC - Boxplot
  output$boxplot <- renderPlotly({
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

  ## QC - Headers - Histogram
  output$counthist <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Count data distributions - histogram")
    }
  })

  # QC - Histogram
  output$hist <- renderPlotly({
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

  ## QC - Headers - Barplot
  output$counttotal <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Total reads")
    }
  })

  # QC - Barplot
  output$barplot <- renderPlotly({
    input$goqc
    tmp <- ddstran()[[1]]
    tmp <- assay(tmp)
    lab <- ddstran()[[2]]
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
  })

  ## QC - Headers - PCA
  output$countpca <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h4("Principal component analysis")
    }
  })

  ## QC - Choose factor - PCA
  output$pcafact <- renderUI({
    # input$goqc
    tmp <- ddsout()[[2]]
    selectInput(
      inputId = "pcafact",
      label = "Choose factor",
      choices = colnames(tmp)
    )
  })

  # QC - PCA
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



  ### ブ レ ー ク  B R E A K  ブ レ ー ク ###



  # DEG - Conditional levels
  output$degfact <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      tmp <- ddsout()[[2]]
      selectInput("fact", "Define factor", colnames(tmp))
    }
  })
  output$deglev1 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {      
      tmp <- ddsout()[[2]]
      selectInput("lev1", "Define level 1", tmp[, input$fact])
    }
  })
  output$deglev2 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      tmp <- ddsout()[[2]]
      selectInput("lev2", "Define level 2", tmp[, input$fact])
    }
  })
  output$deseqform <- renderUI({
    if (input$dgetype != "deseq") {
      return()
    } else {
      textInput(
        inputId = "deseqform",
        label = "Enter experimental design",
        value = "" 
      )
    }
  })

  deg <- eventReactive(input$godeg, {
    if (input$dgetype == "deseq") {
      withProgress(message = "Running DESeq2 analysis...", value = 0, {
        incProgress()
        dds <- DESeqDataSetFromMatrix(
          countData = ddsout()[[3]],
          colData = ddsout()[[2]],
          design = as.formula(input$deseqform)
        )
        dds <- DESeq(dds)
        incProgress()
        fact <- input$fact
        lev1 <- input$lev1
        lev2 <- input$lev2
        incProgress()
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- subset(datobj, baseMean != 0)
        datobj <- datobj %>% tibble::rownames_to_column()
        names(datobj)[1] <- "id"
        incProgress()
      })
    } else if (input$dgetype == "limma") {
      withProgress(message = "Running limma-voom analysis...", value = 0, {
        fact <- input$fact
        lev1 <- input$lev1
        lev2 <- input$lev2
        coldata <- ddsout()[[2]]
        cts <- ddsout()[[3]]
        design <- model.matrix(~ 0 + coldata[, fact])
        rownames(design) <- rownames(coldata)
        colnames(design) <- levels(coldata[, fact])
        cont.des <- paste(lev1, "-", lev2)
        head <- paste0(lev1, "VS", lev2)
        cont <- makeContrasts(contrasts = cont.des, levels = design)
        colnames(cont) <- head
        incProgress()
        v <- voom(cts, design)
        incProgress()
        fit <- lmFit(v)
        incProgress()
        fit.cont <- contrasts.fit(fit, cont)
        fit.cont <- eBayes(fit.cont)
        de.genes <- topTable(fit.cont, coef = 1, number = nrow(cts))
        datobj <- de.genes[order(row.names(de.genes)), ]
        datobj$baseMean <- rowMeans(cts)
        names(datobj) <- c(
          "log2FoldChange",
          "avgexpr",
          "t",
          "pvalue",
          "padj",
          "B",
          "baseMean"
        )
        datobj <- subset(datobj, baseMean != 0)
        datobj <- datobj %>% tibble::rownames_to_column()
        names(datobj)[1] <- "id"
        incProgress()
      })
    } else if (input$dgetype == "edger") {
      withProgress(message = "Running edgeR analysis...", value = 0, {
        incProgress()
        fact <- input$fact
        lev1 <- input$lev1
        lev2 <- input$lev2
        coldata <- ddsout()[[2]]
        cts <- ddsout()[[3]]
        design <- model.matrix(~ 0 + coldata[, fact])
        rownames(design) <- rownames(coldata)
        colnames(design) <- levels(coldata[, fact])
        cont.des <- paste(lev1, "-", lev2)
        head <- paste0(lev1, "VS", lev2)
        cont <- makeContrasts(contrasts = cont.des, levels = design)
        colnames(cont) <- head
        incProgress()
        dge <- DGEList(counts = cts)
        dge <- calcNormFactors(dge, method = "TMM")
        incProgress()
        dge <- estimateGLMCommonDisp(dge, design)
        incProgress()
        dge <- estimateGLMTrendedDisp(dge, design)
        incProgress()
        dge <- estimateGLMTagwiseDisp(dge,design)
        incProgress()
        fit.edger <- glmFit(dge, design)
        lrt.edger <- glmLRT(fit.edger, contrast = cont)
        res <- topTags(lrt.edger, n = nrow(pas.cts), sort.by = "none")
        res <- as.data.frame(res)
        datobj <- res[order(row.names(res)), ]
        datobj$baseMean <- rowMeans(cts)
        names(datobj) <- c(
          "log2FoldChange",
          "logCPM",
          "LR",
          "pvalue",
          "padj",
          "baseMean"
        )
        datobj <- subset(datobj, baseMean != 0)
        datobj <- datobj %>% tibble::rownames_to_column()
        names(datobj)[1] <- "id"
        incProgress()
      })
    }
    return(list(datobj))
  })

  # output$debugdeg <- renderPrint({
  #   debug1 <- deg()[[1]]
  #   head(debug1)
  # })

  ## DEG - Headers - Conditional header
  output$degcomp <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      h2("Comparisons")
    } 
  })

  ## DEG - Headers - Data load check
  output$degcheck <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      p("Data comparisons have been successfully loaded. Please submit visualization parameters.")
    } 
  })


  ## DEG - Headers - Visualization parameters
  output$degvishead <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      h4("Visualization Parameters")
    } 
  })

  ## DEG - select input - Visualization type
  output$vistype <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      selectInput(
        inputId = "plottype",
        label = "Choose plot type",
        choices = c(
          "MA plot" = "maplot",
          "Volcano plot" = "volplot"
        ),
        selected = ""
      )
    }
  })

  ## DEG - Filters - LFC
  output$deglfcfilt <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      textInput(
        inputId = "lfcfilt",
        label = "Choose abs(LFC) cutoff",
        value = 0
      )
    }
  })

  ## DEG - Filters - LFC
  output$degpadfilt <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      textInput(
        inputId = "padfilt",
        label = "Choose padj cutoff",
        value = 1
      )
    }
  })

  ## DEG - Action button - submit vis filters
  output$govisfilt <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      actionButton("govisfilt", "Submit", icon = icon("magic"))
    } 
  })
  
  # DEG - Filter Data 
  deg2 <- eventReactive(input$govisfilt, {
    tmp <- deg()[[1]]
    tmp$padj <- round(tmp$padj, 3)
    tmp[is.na(tmp)] <- 1
    padj <- input$padfilt
    lfc <- input$lfcfilt
    tmp <- tmp[abs(tmp$log2FoldChange) >= lfc, ]
    tmp <- tmp[tmp$padj <= padj, ]
    return(tmp)
  })



  # DEG - Shared data
  share <- eventReactive(input$govisfilt, {
    SharedData$new(deg2())
  })



  # DEG - Conditional plot
  output$degplot <- renderPlotly({
    s <- input$mytable_rows_selected
    datobj <- deg2()
    tooltips <- paste0(
      "<b>ID:</b> ", datobj$id, "<br />",
      "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
      "<b>PADJ:</b> ", round(datobj$padj, 3), "<br />",
      "<b>BM:</b> ", round(datobj$baseMean, 3) 
    )
    if (input$godeg == 0) {
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
          pp <- deg2() %>%
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
            data = deg2()[s, , drop = FALSE], 
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
          pp <- deg2() %>%
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
            data = deg2()[s, , drop = FALSE], 
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
  })


  # DEG - Conditional table
  output$mytable <- DT::renderDataTable({
    if (input$godeg == 0) {
      return()
    } else {
      tmp <- deg2() %>% mutate_if(is.numeric, round, digits = 3)
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
        expr = !is.null(input$govisfilt),
        message = ""   
      )
    )
    if(input$govisfilt == 0) {
      return()
    } else {
      downloadButton("downfiltData", "Download Filtered Data")
    }     
  })

  # Download FILTERED data
  output$downfiltData <- downloadHandler(
    filename = function() {
      paste(
        "filtered-results-", 
        nrow(deg2()[input[["mytable_rows_all"]], ]),
        ".csv", 
        sep = ""
      )
    },
    content = function(file) {
      write.csv(
        deg2()[input[["mytable_rows_all"]], ], 
        file, row.names = FALSE
      )
    }
  )   

  # DEG - Show download button (ALL)
  output$downloadall <- renderUI({
    validate(
      need(
        expr = !is.null(input$govisfilt),
        message = ""   
      )
    )
    if(input$govisfilt == 0) {
      return()
    } else {
      downloadButton("downallData", "Download All Data")
    }     
  })

  # Download ALL data
  output$downallData <- downloadHandler(
    filename = function() {
      paste("all-results.csv")
    },
    content = function(file) {
      write.csv(
        deg()[[1]], 
        file, row.names = FALSE
      )
    }
  )   



  ### ブ レ ー ク  B R E A K  ブ レ ー ク ###



  # Heat - Header
  output$heattitle <- renderUI({
    if(input$goheat == 0) {
      return()
    } else {
      h2(paste0("Heatmap of the ", input$heatrows, " most variable IDs"))
    }
  })

  heattran <- eventReactive(input$goheat, {
    dds <- ddsout()[[1]]
    if (input$heattransform == "log") {
      tmp <- normTransform(dds)
      tmp <- assay(tmp)
      lab <- "log<sub>2</sub>(counts + 1)"
    } else if (input$heattransform == "rlog") {
      tmp <- rlog(dds)
      tmp <- assay(tmp)
      lab <- "rlog(counts)"
    } else if (input$heattransform == "vst") {
      tmp <- vst(dds)
      tmp <- assay(tmp)
      lab <- "vst(counts)"
    } else if (input$heattransform == "raw") {
      tmp <- dds
      tmp <- counts(dds)
      lab <- "Raw counts"
    }
    num <- input$heatrows
    return(list(tmp, lab, num))    
  })

  # heattran2 <- eventReactive(input$goheat, {
  #   heat.counts <- heattran()[[1]]
  #   heat.res <- deg()[[1]]
  #   names(heat.res)[1] <- "row"
  #   heat.res <- as_tibble(heat.res)
  #   heat <- heat.counts[arrange(heat.res, padj, pvalue)$row[1:10], ]
  #   heat <- t(heat)
  #   heat <- scale(heat)
  #   heat <- t(heat)
  #   return(list(heat.counts, heat.res, heat))
  # })





  # Heat - The actual heatmap
  output$heatplot <- renderPlotly({
    heat.counts <- heattran()[[1]]
    heat.res <- deg()[[1]]
    names(heat.res)[1] <- "row"
    heat.res <- as_tibble(heat.res)
    num <- heattran()[[3]]
    heat <- heat.counts[arrange(heat.res, padj, pvalue)$row[1:num], ]
    heat <- t(heat)
    heat <- scale(heat)
    heat <- t(heat)
    tooltips <- paste0(
      "<b>Sample:</b> ", colnames(heat), "<br />",
      "<b>ID:</b> ", rownames(heat), "<br />",
      "<b>Count:</b> ", round(heat, 3)
    )
    tooltips <- matrix(tooltips, ncol = ncol(heat), byrow = TRUE)          
    plot_ly(
      x = colnames(heat),
      y = rownames(heat),
      z = heat,
      type = "heatmap",
      text = tooltips,
      hoverinfo = "text",
      source = "heatplot"
    ) %>%
    layout(
      xaxis = list(title = ""),
      yaxis = list(title = ""),
      margin = list(l = 100)
    )
  })


  output$heatcount <- renderPlotly({
    if (input$goheat == 0) {
      return()
    } else {
      s <- event_data("plotly_click", source = "heatplot")
      source("misc-functions.R")
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
        x = test[, input$fact],
        y = test[, "counts"],
        color = test[, input$fact],
        text = tooltips,
        hoverinfo = "text" 
      ) %>%
      layout(
        title = paste(s[["y"]], "Counts"),
        xaxis = list(title = paste(input$fact)),
        yaxis = list(title = "Normalized counts")
      )
    }
  })
}