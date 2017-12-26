## ViDGER Tutorial

### Accessibility
ViDGER can be freely accessed directly through (**ADD LINK**) or through R using the following commands:

```{r}
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("shiny-tests", "btmonier", subdir = "13-vidger-final")
```

Typically, the link will provide an easier route to using ViDGER.  In circumstances where internet connections will be limited (such as during travel), loading ViDGER through R while internet is still available will allow users to utilize ViDGER without an internet connection later on.

### Input Data
ViDGER requires two pieces of information for analysis.  The first is an expression estimation matrix, also referred to as a count matrix, displaying the gene expression estimates for each sample.  The format requires a CSV file with the row names to list the gene IDs and column names to list the sample IDs.  The second required input is a condition matrix, wherein the factor levels for each sample are provided.  This file requires a CSV format and row names to be the sample IDs matching the sample IDs from the expression estimation matrix and the column names to be the condition factors.

The data used for this tutorial are derived from 28 *Vitis vinifera* samples with three distinct factors (Rootstock, row, and block).

**Expression Matrix (first 4 samples IDs)**

$$\sum_{i=1}^{n}\left( \frac{X_i}{Y_i} \right)$$


