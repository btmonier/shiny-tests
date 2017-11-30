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
      datobj <- lfcShrink(dds, contrast = c(facct, lev1, lev2), res = res)
      datobj <- as.data.frame(datobj)
      datobj <- subset(datobj, baseMean != 0)
      datobj <- datobj %>% tibble::rownames_to_column()
      incProgress()
      return(list(datobj, dds))
    })
  } else if (input$dgetype == "edger") {
    withProgress(message = "Running edgeR analysis...", value = 0, {
      fact <- input$fact
      lev1 <- input$lev1
      lev2 <- input$lev2
      coldata <- ddsout()[[3]]
      cts <- ddsout()[[2]]
      design <- model.matrix(~ 0 + coldata[, fact])
      rownames(design) <- rownames(coldata)
      colnames(design) <- levels(coldata[, fact])
      cont <- paste(lev1, "-", lev2)
      head <- paste0(lev1, "VS", lev2)
      cont <- makeContrasts(cont, levels = design)
      colnames(cont) <- head
      incProgress()
      v <- voom(cts, design)
      incProgress()
      fit <- lmFit(v)
      incProgress()
      fit.cont <- contrasts.fit(fit, cont)
      fit.cont <- eBayes(fit.cont)
      de.genes <- topTable(fit.cont, coef = 1, number = nrow(pas.cts))
      de.genes <- de.genes[order(row.names(de.genes)), ]
      return(list(de.genes))
    })
  } else if (input$dgetyp == "limma") {
    withProgress(message = "Running limma-voom analysis...", value = 0, {

    })
  }
})