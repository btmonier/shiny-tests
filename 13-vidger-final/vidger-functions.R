#------------------------------------------------------------------------------
# Title:  Shiny Test 12 - Experimental Designs - Miscellaneous Functions
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   12.08.17
#------------------------------------------------------------------------------

# HOUSE KEEPING FUNCTIONS

## CRAN Install
pack.man <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)
}

## Size factor extraction (for normalization)
getSizeFact <- function(rc.data) {
	geomMean <- function(x) { 
		prod(x)^(1 / length(x)) 
	}
	gm.mean <- apply(rc.data, 1, geomMean)
	gm.mean[gm.mean == 0] <- NA
	rc.data <- sweep(rc.data, 1, gm.mean, FUN = "/") 
	sf <- apply(rc.data, 2, median, na.rm = TRUE)
	return(sf)
}

## Get normalized counts from different object classes
getNormCounts <- function(rc.data) {
  if (class(rc.data) == "DESeqDataSet") {
    nc.data <- BiocGenerics::counts(rc.data, normalize = TRUE)
    return(nc.data)
  } else {
    sf <- getSizeFact(rc.data)
    nc.data <- sweep(rc.data, 2, sf, FUN = "/")
    return(nc.data)
  }
}

## Get normalized counts for selected gene (heatmap interactivity)
getGenes <- function(rc.data, id, coldata) {
  nc.data <- getNormCounts(rc.data)
  nc.data <- as.data.frame(nc.data[rownames(nc.data) == id, ])
  names(nc.data) <- "counts"
  dat.l <- list(coldata, nc.data)
  nc.data <- Reduce(
    merge, lapply(dat.l, function(x) data.frame(x, sample = row.names(x)))
  )
  return(nc.data)
}

## Get contrast tables from different object classes
getContTable <- function(de.genes, coef, cts, expset, design, fact) {
  if (class(de.genes) == "MArrayLM") {
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
  } else if (class(de.genes) == "DGEGLM") {
    if (expset == "exp1" | expset == "exp2") {
      de.genes2 <- glmLRT(
        glmfit = de.genes,
        contrast = design[, coef]
      )
      de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "logCPM",
        "LR",
        "pvalue",
        "padj",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"
    } else if (expset == "exp3" | expset == "exp4") {
      de.genes2 <- glmLRT(
        glmfit = de.genes,
        coef = coef
      )
      de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2$baseMean <- rowMeans(cts)
      names(de.genes2) <- c(
        "log2FoldChange",
        "logCPM",
        "LR",
        "pvalue",
        "padj",
        "baseMean"
      )
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()      
    }
  } else if (class(de.genes) == "DESeqDataSet") {
    if (expset == "exp1") {
      coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
      coef2 <- unlist(coef2)
      de.genes2 <- results(
        object = de.genes,
        contrast = c(fact, coef2[1], coef2[2])
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"      
    } else if (expset == "exp2") {
      coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
      coef2 <- unlist(coef2)
      de.genes2 <- results(
        object = de.genes,
        contrast = c("group", coef2[1], coef2[2])
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"         
    } else if (expset == "exp3" | expset == "exp4") {
      names(mcols(de.genes))[grep(
        "log2 fold change", 
        mcols(mcols(de.genes))$description
      )] <- colnames(design)
      de.genes2 <- results(
        object = de.genes,
        contrast = list(coef)
      )
      de.genes2 <- as.data.frame(de.genes2)
      de.genes2 <- subset(de.genes2, baseMean != 0)
      de.genes2 <- de.genes2 %>% tibble::rownames_to_column()
      names(de.genes2)[1] <- "id"         
    }
  }
  return(de.genes2)
}



# LIMMA RETURN FIT functions

## LIMMA - EXP 1 - two group comparisons
limma.exp1 <- function(fact, coldata, cts, perm.h) {
  design <- model.matrix(~ 0 + coldata[, fact])
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  v <- voom(cts, design)
  fit <- lmFit(v)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  return(list(fit.cont, cont))  
}

## LIMMA - EXP 2 - multiple factor comparisons
limma.exp2 <- function(fact1, fact2, coldata, cts, perm.h) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  design <- model.matrix(~ 0 + group.c)
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(group.c)
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  v <- voom(cts, design)
  fit <- lmFit(v)
  fit.cont <- contrasts.fit(fit, cont)
  fit.cont <- eBayes(fit.cont)
  return(list(fit.cont, cont)) 
}

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



# EDGER RETURN FIT functions

## EDGER - EXP1 - two group comparisons
edger.exp1 <- function(fact, coldata, cts, perm.h, norm) {
  design <- model.matrix(~ 0 + coldata[, fact])
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont)) 
}

## EDGER - EXP2 - multiple factor comparisons
edger.exp2 <- function(fact1, fact2, coldata, cts, perm.h, norm) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  design <- model.matrix(~ 0 + group.c)
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(group.c)
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, cont)) 
}

## EDGER - EXP3 - classical interactions
edger.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl, 
                       norm) {
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
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, design))
}

## EDGER - EXP4 - added effects blocking and paired
edger.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl,
                       norm) {
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
  # design <- cbind(tmp1, tmp2, tmp3)
  dge <- DGEList(counts = cts)
  dge <- calcNormFactors(dge, method = norm)
  dge <- estimateGLMCommonDisp(dge, design)
  dge <- estimateGLMTrendedDisp(dge, design)
  dge <- estimateGLMTagwiseDisp(dge, design)
  fit.edger <- glmFit(dge, design)
  return(list(fit.edger, design))
}



# DESeq2

## DESeq2 - EXP1 - Two group comparisons
deseq.exp1 <- function(fact, coldata, cts, perm.h) {
  design0 <- model.matrix(~ 0 + coldata[, fact])
  rownames(design0) <- rownames(coldata)
  colnames(design0) <- levels(coldata[, fact])
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont0 <- makeContrasts(contrasts = perm.c, levels = design0)
  colnames(cont0) <- perm.h

  design <- paste("~", fact)
  design <- as.formula(design)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design
  )
  dds <- DESeq(dds)

  return(list(dds, cont0))
}

## DESeq2 - EXP2 - multiple factor comparisons
deseq.exp2 <- function(fact1, fact2, coldata, cts, perm.h) {
  group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))
  design <- model.matrix(~ 0 + group.c)
  rownames(design) <- rownames(coldata)
  colnames(design) <- levels(group.c)
  perm.c <- gsub(pattern = "_VS_", replacement = "-", perm.h)
  cont <- makeContrasts(contrasts = perm.c, levels = design)
  colnames(cont) <- perm.h
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = ~ 1
  )
  dds$group <- group.c
  design(dds) <- ~ group
  dds <- DESeq(dds)
  return(list(dds, cont))
}

## EDGER - EXP3 - classical interactions
deseq.exp3 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
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

  design.dds <- paste("~", fact1, "*", fact2)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds)
  return(list(dds, design))
}

## EDGER - EXP4 - added effects blocking and paired
deseq.exp4 <- function(fact1, fact2, coldata, cts, fact1.rlvl, fact2.rlvl) {
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
  # design <- cbind(tmp1, tmp2, tmp3)
  design.dds <- paste("~", fact1, "+", fact2)
  design.dds <- as.formula(design.dds)
  dds <- DESeqDataSetFromMatrix(
    countData = cts,
    colData = coldata,
    design = design.dds
  )
  dds <- DESeq(dds)
  return(list(dds, design))
}



# PLOT Functions

## Bicluster plots
bicPlot <- function(n, res, cts.var) {
  n <- as.numeric(n)
  par(mar = c(10, 4, 3, 5) + 0.1, yaxt = "n")
  quheatmap(
    x = cts.var,
    bicResult = res,
    number = n, 
    showlabel = TRUE
  )  
}