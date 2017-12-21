# Limma debug

## LIMMA - EXP 3 - classical interactions
limma.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] * coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  tmp1 <- design[, 1:f1n, drop = FALSE]
  tmp2 <- design[, (f1n + 1):dim(design)[2], drop = FALSE]
  drops <- grep("\\:", colnames(tmp2))
  tmp3 <- tmp2[, drops, drop = FALSE]
  tmp2 <- tmp2[, -drops, drop = FALSE]
  colnames(tmp1)[-1] <- paste0(colnames(tmp1)[-1], "_VS_", fact1.rlvl)
  colnames(tmp1)[1] <- gsub("\\(|\\)", "", colnames(tmp1)[1])
  colnames(tmp2) <- paste0(colnames(tmp2), "_VS_", fact2.rlvl)
  design <- cbind(tmp1, tmp2, tmp3)
  v <- voom(cts, design)
  fit <- lmFit(cts, design)
  fit.cont <- eBayes(fit)
  return(list(fit.cont, design))
}

## LIMMA - EXP 4 - added effects blocking and paired
limma.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
  f1n <- length(levels(coldata[, fact1]))
  coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
  coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)
  design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
  rownames(design) <- rownames(coldata)
  colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
  colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))
  colnames(design)[1] <- gsub("\\(|\\)", "", colnames(design)[1])
  colnames(design)[-1:-f1n] <- paste0(
    colnames(design)[-1:-f1n],
    "_VS_", 
    fact2.rlvl
  )
  v <- voom(cts, design)
  fit <- lmFit(cts, design)
  fit.cont <- eBayes(fit)
  return(list(fit.cont, design))
}

de.genes <- tmp[[1]]
design <- tmp[[2]]
coef <- colnames(design)[2]


de.genes2 <- topTable(
  fit = de.genes,
  coef = coef,
  number = nrow(cts)
)
de.genes2 <- de.genes2[order(rownames(de.genes2)), ]
de.genes2 <- as.data.frame(de.genes2)
de.genes2$baseMean <- rowMeans(cts)
names(de.genes2) <- c(
  "log2FoldChange",
  "avgexpr",
  "t",
  "pvalue",
  "padj",
  "B",
  "baseMean"
)
de.genes2 <- subset(de.genes2, baseMean != 0)
de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
names(de.genes2)[1] <- "id"