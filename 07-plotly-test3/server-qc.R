#------------------------------------------------------------------------------
# Title:  Shiny Test 07 - Modularity - Server Logic (QC)
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.26.17
#------------------------------------------------------------------------------

qc <- function(input, output) {
	reactive({
    if (input$go == 0) {
      return()
    } else {
      isolate({
        cts <- input$file1
        coldata <- input$file2
        cts <- as.matrix(read.csv(cts$datapath, header = TRUE, row.names = 1))
        coldata <- read.csv(coldata$datapath, header = TRUE, row.names = 1)
        cts <- cts[, rownames(coldata)]
        
        if (input$expdes == "") {
          expdes <- NULL
        } else {
          expdes <- input$expdes
        }
        
        dds <- DESeqDataSetFromMatrix(countData = cts,
                                      colData = coldata,
                                      design = as.formula(expdes))
        dds <- dds[ rowSums(counts(dds)) > input$prefilt, ]
        # dds <- DESeq(dds)
        
        # reslfc <- lfcShrink(dds = dds, coef = 2, res = results(dds))
        
        return(list(dds))  
      })
    }
  })
}

qcopttran <- function(input, output) {
  output$opttran <- renderUI({
    validate(
      need(input$go != 0, "")
    )
    selectInput(
      inputId = "transform",
      label = HTML(paste("<br/>Choose transformation method for counts")),
      choices = c(
        "Normal log: log2(n + psuedocount)" = "log",
        "Regularized log: rlog(n)" = "rlog",
        "Variance stabilizing transform: vst(n)" = "vst",
        "No transformation" = "raw"
      )       
    )
  })
}

