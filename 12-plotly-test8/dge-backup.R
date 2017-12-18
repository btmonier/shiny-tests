
  ## DEG - analysis - reactive
  deg <- eventReactive(input$godge, {
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
        res <- topTags(lrt.edger, n = nrow(cts), sort.by = "none")
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