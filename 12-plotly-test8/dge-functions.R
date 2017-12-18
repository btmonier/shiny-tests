  ## DEG - analysis - reactive
  deg <- eventReactive(input$godge, {
    cts <- ddsout()[[3]]
    coldata <- ddsout()[[2]]
    if (input$dgemetho == "limma") {
      if(input$dgeexpsetup == "exp1") {
        fact <- input$dgeexp1a
        design <- model.matrix(~ 0 + coldata[, fact])
        rownames(design) <- rownames(coldata)
        colnames(design) <- levels(coldata[, fact])
        perm.h <- input$dgeexp1b
        perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
        cont <- makeContrasts(contrasts = perm.c, levels = design)
        colnames(cont) <- perm.h
        v <- voom(cts, design)
        fit <- lmFit(v)
        fit.cont <- contrasts.fit(fit, cont)
        fit.cont <- eBayes(fit.cont)
        comp <- perm.h[input$dgemaincontrasts]
        de.genes <- topTable(fit.cont, coef = comp, number = nrow(cts))
        de.genes <- de.genes[order(row.names(de.genes)), ]
      }
    } 
    return(list(de.genes))
  })