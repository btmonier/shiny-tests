#------------------------------------------------------------------------------
# Title:  limma-voom and edgeR Tests
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.29.17
#------------------------------------------------------------------------------

# Preamble ----

## Set working directory
setwd("D:/Box Sync/misc-shiny-apps/00-misc-tests")

## Load data
pas.cts <- as.matrix(read.csv("count-data.csv", header = TRUE, row.names = 1))
pas.coldata <- read.csv("col-data.csv", header = TRUE, row.names = 1)
pas.cts <- pas.cts[, rownames(pas.coldata)]

## Load packages
library(limma)



# Matrices... ----

## Create design and contras matrix (simple models only)
fact <- "condition"
lev1 <- "treated"
lev2 <- "untreated"

design <- model.matrix(~ 0 + pas.coldata[, fact])
rownames(design) <- rownames(pas.coldata)
colnames(design) <- levels(pas.coldata$condition)

cont <- paste(lev1, "-", lev2)
head <- paste0(lev1, "VS", lev2)
cont <- makeContrasts(cont, levels = design)
colnames(cont) <- head



# limma-voom pipeline ----

## Voom transform data
v <- voom(pas.cts, design)


## lmFit and eBayes
fit <- lmFit(v)
fit.cont <- contrasts.fit(fit, cont)
fit.cont <- eBayes(fit.cont)


de.genes <- topTable(fit.cont, coef = 1, number = nrow(pas.cts))
de.genes <- de.genes[order(row.names(de.genes)), ]



# edgeR pipeline ----

## Make DGEList class
dge <- DGEList(counts = pas.cts)
dge <- calcNormFactors(dge, method = "TMM")

## Estimate dispersion parameters for GLM
dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge,design)

## Model fitting
fit.edger <- glmFit(dge, design)

## Differential expression
lrt.edger <- glmLRT(fit.edger, contrast = cont)

## Get results
res <- topTags(lrt.edger, n = nrow(pas.cts), sort.by = "none")
res <- as.data.frame(res)

edger.results <- lrt.edger$table
sig.edger <- decideTestsDGE(lrt.edger, adjust.method = "BH", p.value = 0.05)
genes.edger <- row.names(edger.results)[which(sig.edger != 0)]
