## Welcome to the ViDGER Pipeline

### About
ViDGER (**Vi**isualization of **D**ifferential **G**ene **E**xpression using **R**), is a web-based tool for the analysis of RNA-seq count data. ViDGER is derived from the R package, ViDGER, which can produce information-rich visualizations for the interpretation of differential gene expression (DGE) results from three widely-used tools, *Cuffdiff*, *DESeq2*, and *edgeR*. ViDGER is a **user-friendly** and **interactive** Shiny app for gene expression analysis. This app takes advantage of several popular DGE tools (*DESeq2*, *edgeR*, and *limma*) available through Bioconductor in conjunction with the Plotly and DataTable API libraries for R.

### Submit and QC
To use this app, all you need to upload is two CSV files: a file of raw read counts and another file for sample treatment data. For further information of how these files should look like, take a look at the demonstration data or go to the help section for a detailed walk through.

### Preliminary Analysis
After you have submitted your data, you can analyze the correlation between any treatment sample combination, discover which genes have the highest average expression, plot individual genes and their respective counts, and potentially discover biclusters in your experimental data.

### DGE Analysis
To determine DGE, you currently have the option of using 3 Bioconductor packages mentioned in the prior sections. DGE can be visualized using common techniques, including volcano plots, MA plots, and heat maps. Additionally, DGE data sets can be downloaded as CSV files for further analysis in this section.

### More
Under this tab, you will find information regarding a step-by-step tutorial on how to analyze RNA-seq data using ViDGER, a frequently asked questions page going over common analytical topics, an about us page detailing the creators of this tool and also session info regarding what is used to make this web-application possible.
