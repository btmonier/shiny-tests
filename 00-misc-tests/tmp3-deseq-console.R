#------------------------------------------------------------------------------
# Title:  DESeq2 Console Backchannel
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.06.17
#------------------------------------------------------------------------------

# DESeq2 "standard" pipeline ----

## Load packages
require(DESeq2)
require(ggplot2)
require(airway)
require(plotly)
require(scatterD3)

## Load data
### Pasilla
setwd("D:/Box Sync/misc-shiny-apps/00-misc-tests")

### Airway
data(airway)


## Prep data
### Pasilla
pas.cts <- as.matrix(read.csv("count-data.csv", header = TRUE, row.names = 1))
pas.coldata <- read.csv("col-data.csv", header = TRUE, row.names = 1)
pas.cts <- pas.cts[, rownames(pas.coldata)]
### Airway
air.cts <- assay(airway)
air.coldata <- colData(airway)
air.cts <- air.cts[, rownames(air.coldata)]


## Create DESeq data set
### Pasilla
pas.dds <- DESeqDataSetFromMatrix(
  countData = pas.cts,
  colData = pas.coldata,
  design = ~ type + condition
)
### Airway
air.dds <- DESeqDataSetFromMatrix(
  countData = air.cts,
  colData = air.coldata,
  design = ~ cell
)


## Pre-filtering
pas.dds <- pas.dds[rowSums(counts(pas.dds)) > 1, ]
air.dds <- air.dds[rowSums(counts(air.dds)) > 1, ]


## Run DESeq2
pas.dds <- DESeq(pas.dds)
air.dds <- DESeq(air.dds)


## Results
pas.res <- results(pas.dds)
air.res <- results(air.dds)


## log2 shrink?
pas.reslfc <- lfcShrink(dds = pas.dds, coef = 2, res = results(pas.dds))
air.reslfc <- lfcShrink(dds = air.dds, coef = 2, res = results(air.dds))



# Parameter entry and "priming" ----

## "User" data (USER WILL CHANGE THIS)
dds <- pas.dds
fact <- "condition"
lev1 <- "untreated"
lev2 <- "treated"
padj <- 0.05
lfc <- 2

## Automatic data priming
res <- results(dds, contrast = c(fact, lev1, lev2))
datobj <- lfcShrink(dds, contrast = c(fact, lev1, lev2), res = res)
datobj <- as.data.frame(datobj)
datobj <- subset(datobj, baseMean != 0)
datobj$isDE <- ifelse(datobj$padj < padj, TRUE, FALSE)
datobj$isDE[is.na(datobj$isDE)] <- FALSE
datobj$isLFC <- ifelse(abs(datobj$log2FoldChange) > lfc, TRUE, FALSE)
datobj$isLFC[is.na(datobj$isLFC)] <- FALSE

## Hover metadata
tooltips <- paste0(
  "<b>ID:</b> ", rownames(datobj), "<br />",
  "<b>LFC:</b> ", round(datobj$log2FoldChange, 3), "<br />",
  "<b>FDR:</b> ", round(datobj$padj, 3), "<br />",
  "<b>BM:</b> ", round(datobj$baseMean, 3) 
)


## MA plot (plotly)
plot_ly(
  type = "scatter",
  mode = "markers",
  x = ~log10(datobj$baseMean),
  y = ~datobj$log2FoldChange,
  color = ~datobj$isDE,
  colors = c("#CCCCCC", "#4286f4"),
  text = tooltips,
  hoverinfo = "text"
) %>%
layout(
  xaxis = list(title = "Base mean"),
  yaxis = list(title = "log<sub>2</sub> fold change")
)


## MA plot (scatterD3)
scatterD3(
  x = log10(datobj$baseMean),
  y = datobj$log2FoldChange,
  point_size = 15,
  hover_size = 10,
  col_var = ifelse(datobj$isDE == TRUE, "Yes", "No"),
  colors = c("Yes" = "#4286f4", "No" = "#CCCCCC"),
  symbol_var = ifelse(datobj$isLFC == TRUE, "Yes", "No"),
  tooltip_text = tooltips,
  xlab = "Base mean",
  ylab = "log2 fold change",
  col_lab = paste0("FDR is significant at ", padj),
  symbol_lab = paste0("LFC is larger than ", lfc)
)


## Volcano plot (plotly)
plot_ly(
  type = "scatter",
  mode = "markers",
  x = datobj$log2FoldChange,
  y = -log10(datobj$pvalue),
  color = ~datobj$isDE,
  colors = c("#CCCCCC", "#4286f4"),
  text = tooltips,
  hoverinfo = "text"
) %>%
layout(
  xaxis = list(title = "log<sub>2</sub> fold change"),
  yaxis = list(title = "-log<sub>10</sub>(p-value)")
)


## Volcano plot (scatterD3)
scatterD3(
  x = datobj$log2FoldChange,
  y = -log10(datobj$pvalue),
  point_size = 15,
  hover_size = 10,
  col_var = ifelse(datobj$isDE == TRUE, "Yes", "No"),
  colors = c("Yes" = "#4286f4", "No" = "#CCCCCC"),
  symbol_var = ifelse(datobj$isLFC == TRUE, "Yes", "No"),
  tooltip_text = tooltips,
  xlab = "log2 fold change",
  ylab = "-log10(p-value)",
  col_lab = paste0("FDR is significant at ", padj),
  symbol_lab = paste0("LFC is larger than ", lfc)
)


## Scatter plot

### Transformation conditionals...
col <- as.data.frame(colData(dds))
tran <- "rlog" 
if (tran == "rlog") {
  tmp <- rlog(dds, blind = FALSE)
  tmp <- assay(tmp)
  lab <- "rlog(counts)"
} else if (tran == "vst") {
  tmp <- vst(dds, blind = FALSE)
  tmp <- assay(tmp)
  lab <- "vst(counts)"
} else if (tran == "counts") {
  tmp <- assay(dds)
  lab <- "raw counts"
} else if (tran == "log") {
  tmp <- normTransform(dds)
  tmp <- assay(tmp)
  lab <- "log2(counts + 1)"
}

### Make temp data frame for plotting
nam_x <- row.names(col[which(col[fact] == lev1), ])
nam_y <- row.names(col[which(col[fact] == lev2), ])
x <- rowMeans(tmp[, nam_x])
y <- rowMeans(tmp[, nam_y])

tmp <- data.frame(x, y)

### The plot
scatterD3(
  x = tmp$x,
  y = tmp$y,
  point_size = 15,
  hover_size = 10,
  colors = "#4286f4",
  tooltip_text = tooltips,
  xlab = paste(lev1, lab),
  ylab = paste(lev2, lab)
)



## Heat map (plotly)
num <- 200

### Make data frame
select <- order(
  rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE
)[1:num]

heat <- as.data.frame(tmp)[select, ]
heat <- t(heat)


tooltips <- paste0(
  "<b>ID:</b> ", unlist(colnames(heat)), "<br />",
  "<b>Count:</b> ", t(round(heat, 3)), "<br />",
  "<b>Tran:</b> ", lab, "<br />"
)
tooltips <- matrix(tooltips, ncol = num, byrow = TRUE)

ax <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
)

plot_ly(
  x = 1:num,
  y = rownames(heat),
  z = heat,
  type = "heatmap",
  text = tooltips,
  hoverinfo = "text",
  colorscale = "RdYlBu"
) %>%
layout(
  xaxis = ax
)


## Heat map (d3heatmap)
num <- 50

### Make data frame
select <- order(
  rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE
)[1:num]

heat <- as.data.frame(tmp)[select, ]
heat <- t(heat)

d3heatmap(
  heat,
  scale = "row",
  dendrogram = "none",
  color = "RdYlBu",
  Colv = FALSE
)



## Box plots (plotly)
### Averages for each factor level
source("tmp5-vidger-functions.R")
box <- getDeseqBox(dds, "type", "log")

plot_ly(
  box,
  y = ~value,
  color = ~key,
  colors = "Paired",
  type = "box"
)


### Plot each sample replicate
source("tmp5-vidger-functions.R")
box <- getDeseqBox(dds, "type", "log")

plot_ly(
  box,
  type = "box",
  y = log2(box$value),
  color = ~key
)

tmp <- as.data.frame(assay(dds))
tmp <- tidyr::gather(tmp)

plot_ly(
  tmp,
  type = "box",
  y = log2(tmp$value + 1),
  color = ~key
)


## Remember where you came from
# ia <- read.csv("ia.csv", header = TRUE)
# ggplot(data = ia, aes(x = lng, y = lat, group = 1)) +
#   geom_point()

# plot_ly(
#   data = ia,
#   type = "scatter",
#   mode = "markers",
#   x = ~lng,
#   y = ~lat
# )



## Plot counts (plotly)
gene <- "FBgn0003360" # change

source("tmp5-vidger-functions.R")
(test <- getDeseqGeneCount(dds, "condition", "rlog", gene))

plot_ly(
  data = test,
  type = "scatter",
  mode = "markers",
  x = test[[1]],
  y = ~value,
  color = test[[1]]
)



plot_ly(
  tmp,
  y = ~value,
  color = ~key,
  colors = "Paired",
  type = "box"
)


plot_ly(
  type = "scatter",
  mode = "markers",
  x = tmp4$PC1,
  y = tmp4$PC2,
  color = tmp4$condition,
  symbol = tmp4$type
)