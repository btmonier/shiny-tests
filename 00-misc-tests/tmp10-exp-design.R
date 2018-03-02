#------------------------------------------------------------------------------
# Title:  Experimental Design Workflow (limma and edgeR)
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   12.07.17
#------------------------------------------------------------------------------

# Preamble ----

## Set working directory
setwd("D:/Box Sync/misc-github-repos/shiny-tests/00-misc-tests")

## Load data
pas.cts <- as.matrix(read.csv("count-data.csv", header = TRUE, row.names = 1))
pas.coldata <- read.csv("col-data.csv", header = TRUE, row.names = 1)
pas.cts <- pas.cts[, rownames(pas.coldata)]

## Load packages
library(DESeq2)
library(limma)
library(edgeR)
library(gtools)



# Experimental Design

## Make comparisons of ONE GROUP

### USER INPUTS
cts <- pas.cts
coldata <- pas.coldata
fact <- "condition"

### Create design and contrast matrices
design <- model.matrix(~ 0 + coldata[, fact])
rownames(design) <- rownames(coldata)
colnames(design) <- levels(coldata[, fact])

perm <- levels(coldata[, fact])
perm <- permutations(n = length(perm), r = 2, v = perm)
perm.c <- apply(perm[, 1:2], 1, paste, collapse = "-" )
perm.h <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")

cont <- makeContrasts(contrasts = perm.c, levels = design)
colnames(cont) <- perm.h

### Voom transform data
v <- voom(cts, design)

### lmFit and eBayes
fit <- lmFit(v)
fit.cont <- contrasts.fit(fit, cont)
fit.cont <- eBayes(fit.cont)

### Choose gene tables
comp <- perm.h[2]
de.genes <- topTable(fit.cont, coef = comp, number = nrow(pas.cts))
de.genes <- de.genes[order(row.names(de.genes)), ]



## Make group comparisons of EACH TREATMENT COMBINATION

### USER INPUTS
cts <- pas.cts
coldata <- pas.coldata
fact1 <- "condition"
fact2 <- "type"
group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))

### Create design and contrast matrices
design <- model.matrix(~ 0 + group.c)
rownames(design) <- rownames(coldata)
colnames(design) <- levels(group.c)

perm <- levels(group.c)
perm <- permutations(n = length(perm), r = 2, v = perm)
perm.c <- apply(perm[, 1:2], 1, paste, collapse = "-")
perm.h <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")

cont <- makeContrasts(contrasts = perm.c, levels = design)
colnames(cont) <- perm.h

### Voom transform data
v <- voom(cts, design)

### lmFit and eBayes
fit <- lmFit(v)
fit.cont <- contrasts.fit(fit, cont)
fit.cont <- eBayes(fit.cont)

### Choose gene tables
comp <- perm.h[2]
de.genes <- topTable(fit.cont, coef = comp, number = nrow(pas.cts))
de.genes <- de.genes[order(row.names(de.genes)), ]



## Establish design with INTERACTIONS!

### USER INPUTS
cts <- pas.cts
cts <- cts[rowSums(cts) > 10, ]
coldata <- pas.coldata
fact1 <- "condition"
fact2 <- "type"
fact1.rlvl <- "untreated" 
fact2.rlvl <- "single.read"
f1n <- length(levels(coldata[, fact1]))

### Create baseline references for model
coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)

### Create design and contrast matrices
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

### Voom transform data
v <- voom(cts, design)

### lmFit and eBayes
fit <- lmFit(cts, design)
fit.cont <- eBayes(fit)

### Choose gene tables
comp <- colnames(design)[2]
de.genes <- topTable(fit.cont, coef = comp, number = nrow(pas.cts), adjust.method = "bonferroni")
de.genes <- de.genes[order(row.names(de.genes)), ]

de.genes.f <- de.genes[de.genes$logFC >= 1, ]
de.genes.f <- de.genes[de.genes$adj.P.Val <= 0.05, ]

## Designs with PAIRING OR BLOCKING

### USER INPUTS
cts <- pas.cts
coldata <- pas.coldata
fact1 <- "type" # CHOOSE BLOCK, SAMPLE, INDIVIDUAL, etc.
fact2 <- "condition" # CHOOSE TREATMENT
fact1.rlvl <- "single.read" # BASE REFERENCE (FACTOR 1)
fact2.rlvl <- "untreated" # BASE REFERENCE (FACTOR 2)
f1n <- length(levels(coldata[, fact1]))

### Create baseline references for model
coldata[, fact1] <- relevel(x = coldata[, fact1], ref = fact1.rlvl)
coldata[, fact2] <- relevel(x = coldata[, fact2], ref = fact2.rlvl)

### Create design and contrast matrices
design <- model.matrix(~ coldata[, fact1] + coldata[, fact2])
rownames(design) <- rownames(coldata)
colnames(design) <- gsub("coldata\\[, fact1\\]*", "", colnames(design))
colnames(design) <- gsub("coldata\\[, fact2\\]*", "", colnames(design))

colnames(design)[-1:-f1n] <- paste0(
	colnames(design)[-1:-f1n],
	"_VS_", 
	fact2.rlvl
)

### Voom transform data
v <- voom(cts, design)

### lmFit and eBayes
fit <- lmFit(cts, design)
fit.cont <- eBayes(fit)

### Choose gene tables
comp <- colnames(design)[3]
de.genes <- topTable(fit.cont, coef = comp, number = nrow(pas.cts))
de.genes <- de.genes[order(row.names(de.genes)), ]






# EDGER experiments - interactions
cts <- pas.cts
coldata <- pas.coldata
fact1 <- "condition"
fact2 <- "type"
fact1.rlvl <- "untreated" 
fact2.rlvl <- "single.read"

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
lrt.edger <- glmLRT(fit.edger, contrast = design[, 1])
fit.names <- NULL


tmp <- edger.exp3(
	cts = pas.cts,
	coldata = pas.coldata,
	fact1 = "condition",
	fact2 = "type",
	fact1.rlvl = "untreated",
	fact2.rlvl = "single.read"
)



## Make group comparisons of EACH TREATMENT COMBINATION

### USER INPUTS
cts <- pas.cts
coldata <- pas.coldata
fact1 <- "condition"
fact2 <- "type"
group.c <- factor(paste(coldata[, fact1], coldata[, fact2], sep = "_"))

### Create design and contrast matrices
design <- model.matrix(~ 0 + group.c)
rownames(design) <- rownames(coldata)
colnames(design) <- levels(group.c)

perm <- levels(group.c)
perm <- permutations(n = length(perm), r = 2, v = perm)
perm.c <- apply(perm[, 1:2], 1, paste, collapse = "-")
perm.h <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")

cont <- makeContrasts(contrasts = perm.c, levels = design)
colnames(cont) <- perm.h

coef <- colnames(cont)[1]


dge <- DGEList(counts = cts)
dge <- calcNormFactors(dge, method = norm)
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)
fit.edger <- glmFit(dge, design)
lrt.edger <- glmLRT(fit.edger, coef = coef)
fit.names <- NULL


test.fun <- function(de.genes, coef, cts, expset, fact) {
	if (class(de.genes) == "MArrayLM") {
    de.genes2 <- topTable(
      fit = de.genes,
      coef = coef,
      number = nrow(cts)
    )
    de.genes2 <- de.genes2[order(rownames(de.genes2)), ]
	} else if (class(de.genes) == "DGEGLM") {
		if (expset == "exp1" | expset == "exp2") {
			de.genes2 <- glmLRT(
				glmfit = de.genes,
				contrast = coef
			)
	    de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
		} else if (expset == "exp3" | expset == "exp4") {
			de.genes2 <- glmLRT(
				glmfit = de.genes,
				coef = coef
			)
			de.genes2 <- topTags(de.genes2, n = nrow(cts), sort.by = "none")
      de.genes2 <- as.data.frame(de.genes2)
		}
	} else if (class(de.genes) == "DESeqDataSet") {
		coef2 <- strsplit(coef, "_VS_", fixed = TRUE)
		coef2 <- unlist(coef2)
		de.genes2 <- results(
			object = de.genes,
			contrast = c(fact, coef2[1], coef2[2])
		)
		de.genes2 <- as.data.frame(de.genes2)
	}
	return(de.genes2)
}



# DESeq2 exp design layout

## Two group comparisons

### USER INPUTS
cts <- pas.cts
coldata <- pas.coldata
fact <- "condition"

### Create design and contrast matrices
design <- model.matrix(~ 0 + coldata[, fact])
rownames(design) <- rownames(coldata)
colnames(design) <- levels(coldata[, fact])

perm <- levels(coldata[, fact])
perm <- permutations(n = length(perm), r = 2, v = perm)
perm.c <- apply(perm[, 1:2], 1, paste, collapse = "-" )
perm.h <- apply(perm[, 1:2], 1, paste, collapse = "_VS_")

cont <- makeContrasts(contrasts = perm.c, levels = design)
colnames(cont) <- perm.h

perm.h <- c("treated_VS_untreated", "untreated_VS_treated")
deseq.exp1("condition", pas.coldata, pas.cts, perm.h)


names(mcols(dds))[grep("log2 fold change", mcols(mcols(dds))$description)] <-
colnames(design)


deseq.exp3(fact1 = "condition", fact2 = "type", coldata = coldata, cts = cts, fact1.rlvl = "untreated", fact2.rlvl = "single.read")