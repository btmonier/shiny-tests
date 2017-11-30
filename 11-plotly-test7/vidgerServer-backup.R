#------------------------------------------------------------------------------
# Title:  Shiny Test 11 - Crosstalk - Server Logic
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.24.17
#------------------------------------------------------------------------------

vidgerServer <- function(input, output) {

	# QC - Reactive
	ddsout <- eventReactive(input$goqc, {
    cts <- input$file1
    coldata <- input$file2

    cts <- as.matrix(read.csv(cts$datapath, header = TRUE, row.names = 1))
    coldata <- read.csv(coldata$datapath, header = TRUE, row.names = 1)
    cts <- cts[, rownames(coldata)]
    
    dds <- DESeqDataSetFromMatrix(
      countData = cts,
      colData = coldata,
      design = ~ 1
    )
    dds <- dds[ rowSums(counts(dds)) > input$prefilt, ]
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

  # QC - Header
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
    plot_ly(
      data = pca,
      type = "scatter",
      mode = "markers",
      x = ~PC1,
      y = ~PC2,
      color = ~group,
      marker = list(size = 8)
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
      tmp <- input$file2
      tmp <- read.csv(tmp$datapath, header = TRUE, row.names = 1)
      selectInput("fact", "Define factor", colnames(tmp))
    }
  })
  output$deglev1 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {      
      tmp <- input$file2
      tmp <- read.csv(tmp$datapath, header = TRUE, row.names = 1)
      selectInput("lev1", "Define level 1", tmp[, input$fact])
    }
  })
  output$deglev2 <- renderUI({
    if (input$goqc == 0) {
      return()
    } else {
      tmp <- input$file2
      tmp <- read.csv(tmp$datapath, header = TRUE, row.names = 1)
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

  # Reactive DEG object
  deg <- reactive({
    if (input$godeg == 0) {
      return()
    } else {
      withProgress(
        message = "Running DEG analysis", 
        detail = "<br/>This may take several minutes...", 
        value = 0, {
        incProgress()
        dds <- ddsout()[[1]]
        dds <- DESeq(dds)
        ## Establish objects
        incProgress()
        fact <- input$fact
        lev1 <- input$lev1
        lev2 <- input$lev2
        ## Automatic data priming
        incProgress()
        res <- results(dds, contrast = c(fact, lev1, lev2))
        datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
        datobj <- as.data.frame(datobj)
        datobj <- subset(datobj, baseMean != 0)
        datobj <- datobj %>% tibble::rownames_to_column()
        incProgress()
        return(list(datobj, dds))
      })
    }
  })

  ## DEG - Headers - Conditional header
  output$degcomp <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
    	h4("Comparisons")
    } 
  })

  # DEG - Filter Data 
  deg2 <- reactive({
    if (input$godeg == 0) {
      return()
    } else {
      tmp <- deg()[[1]]
      # tmp <- tmp %>% mutate_if(is.numeric, round, digits = 3)
      tmp$padj <- round(tmp$padj, 3)
      tmp[is.na(tmp)] <- 1
      padj <- input$padfilt
      lfc <- input$lfcfilt
      tmp <- tmp[abs(tmp$log2FoldChange) >= lfc, ]
      tmp <- tmp[tmp$padj <= padj, ]
      return(tmp)
    }
  })

  # DEG - Shared data
  share <- reactive({
    SharedData$new(deg2())
  })

  ## DEG - Headers - Visualization parameters
  output$vishead <- renderUI({
    if(input$godeg == 0) {
      return()
    } else {
      h4("Visualization Parameters")
    }
  })

  ## DEG - Headers - Visualization type
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



  # DEG - Conditional plot
  output$degplot <- renderPlotly({
    s <- input$mytable_rows_selected
    if (input$godeg == 0) {
      return()
    } else {
    	# isolate({
    		if (input$plottype == "maplot") {
    		  if (!length(s)) {
    		    p <- share() %>%
    		      plot_ly(
    		      	x = ~log10(baseMean), 
    		      	y = ~log2FoldChange,
                type = "scatter", 
    		      	mode = "markers", 
    		      	color = I("darkgray"), 
    		      	name = "Unfiltered") %>%
    		      layout(showlegend = TRUE) %>% 
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
    		      	name = "Unfiltered") %>%
    		      layout(showlegend = TRUE)
		
    		    # selected data
    		    pp <- add_trace(
    		    	pp, 
    		    	data = deg2()[s, , drop = FALSE], 
    		    	x = ~log10(baseMean), 
    		    	y = ~log2FoldChange,
              type = "scatter", 
    		    	mode = "markers", 
    		    	color = I("royalblue1"), 
    		    	size = I(10), 
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
    		      	name = "Unfiltered") %>%
    		      layout(showlegend = TRUE) %>% 
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
    		      	name = "Unfiltered") %>%
    		      layout(showlegend = TRUE)
		
    		    # selected data
    		    pp <- add_trace(
    		    	pp, 
    		    	data = deg2()[s, , drop = FALSE], 
    		    	x = ~log2FoldChange, 
    		    	y = ~-log10(pvalue), 
    		    	mode = "markers", 
    		    	color = I("royalblue1"), 
    		    	size = I(10), 
    		    	name = "Filtered"
    		    )
    		  }      
    		}    		
    	# })
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
      	"rowname", 
      	target = "row", 
      	color = DT::styleEqual(m2$rowname, rep("white", length(m2$rowname))),
      	backgroundColor = DT::styleEqual(
      		m2$rowname, 
      		rep("darkgray", length(m2$rowname))
      	)
      )
    }
  }
  })  



  ### B R E A K ###



  # Heat - Header
  output$heattitle <- renderUI({
    if(input$goqc == 0) {
      return()
    } else {
      h2(paste0("Heatmap of the ", input$heatrows, " most DEGs"))
    }
  })

  # Heat - The actual heatmap
  output$heatplot <- renderPlotly({
    if(input$goqc == 0 && input$godeg == 0) {
      return()
    } else {
      if (input$transform == "log") {
        tmp <- deg()[[2]]
        tmp <- normTransform(dds)
        tmp <- assay(tmp)
        lab <- "log<sub>2</sub>(counts + 1)"
      } else if (input$transform == "rlog") {
        tmp <- deg()[[2]]
        tmp <- rlog(dds)
        tmp <- assay(tmp)
        lab <- "rlog(counts)"
      } else if (input$transform == "vst") {
        tmp <- deg()[[2]]
        tmp <- vst(dds)
        tmp <- assay(tmp)
        lab <- "vst(counts)"
      } else if (input$transform == "raw") {
        tmp <- deg()[[2]]
        tmp <- assay(dds)
        lab <- "Raw counts"
      }

      dds <- deg()[[2]]
      res <- results(
        dds, 
        contrast = c(input$fact, input$lev1, input$lev2), 
        tidy = TRUE
      )
      res <- as_tibble(res)
      
      num <- input$heatrows
      heat <- tmp[arrange(res, padj, pvalue)$row[1:num], ]
      heat <- t(heat)
      heat <- scale(heat)
      heat <- t(heat)      
      plot_ly(
        x = colnames(heat),
        y = rownames(heat),
        z = heat,
        type = "heatmap",
        source = "heatplot"
      ) %>%
      layout(
        xaxis = list(title = ""),
        yaxis = list(title = ""),
        margin = list(l = 100)
      )
    }
  })


  output$heatcount <- renderPlotly({
    s <- event_data("plotly_click", source = "heatplot")
    source("misc-functions.R")
    test <- getGenes(
      rc.data = deg()[[2]], 
      id = s[["y"]],
      coldata = ddsout()[[2]]
    )
    plot_ly(
      data = test,
      type = "scatter",
      mode = "markers",
      x = test[, input$fact],
      y = test[, "counts"],
      color = test[, input$fact]
    ) %>%
    layout(
      title = paste(s[["y"]], "Counts"),
      xaxis = list(title = paste(input$fact)),
      yaxis = list(title = "Normalized counts")
    )   
  })
}