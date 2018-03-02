#------------------------------------------------------------------------------
# Title:  User Input Tests for Model Matrices
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   2018-02-28
#------------------------------------------------------------------------------

# Load data
setwd("D:/Box Sync/misc-github-repos/shiny-tests/00-misc-tests")
pas.cts <- as.matrix(read.csv("count-data.csv", header = TRUE, row.names = 1))
pas.coldata <- read.csv("col-data.csv", header = TRUE, row.names = 1)
pas.cts <- pas.cts[, rownames(pas.coldata)]


# Model matrix (not)
condition <- pas.coldata[, "condition"]
type <- pas.coldata[, "type"]
design <- model.matrix(~ condition * type)
rownames(design) <- rownames(pas.coldata)

# Pass model matrix into DESeq2...
pas.dds <- DESeqDataSetFromMatrix(
  countData = pas.cts,
  colData = pas.coldata,
  design = ~ 1
)
pas.dds <- DESeq(pas.dds, full = design, modelMatrixType = "standard")


# With limma-voom...

fit <- lmFit(pas.cts, design)
fit <- eBayes(fit)
de.genes2 <- topTable(
  fit = fit,
  coef = 3,
  number = nrow(pas.cts)
)
de.genes2 <- de.genes2[order(rownames(de.genes2)), ]
de.genes2 <- as.data.frame(de.genes2)
de.genes2$baseMean <- rowMeans(pas.cts)
names(de.genes2) <- c(
  "log2FoldChange",
  "avgexpr",
  "t",
  "pvalue",
  "padj",
  "B",
  "baseMean"
)

de.genes2[de.genes2$padj <= 0.05, ]