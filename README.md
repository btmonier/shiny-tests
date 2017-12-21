# shiny-tests

## Overview

This is a collection of Shiny tests to determine the interactivity of the `vidger` package.

To run the examples locally, you can install the shiny package in R, and use the function `runGithub()`. For example, to run `13-vidger-final`:

``` r
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("shiny-tests", "btmonier", subdir = "13-vidger-final")
```


## Fix list

```
Last updated: 2017-12-20 21:34:02 CST
```

| Task                                | Completed?     |
|-------------------------------------|----------------|
| Add clustering, heatmap to "Submit" | X					     |
| Add DEG Analysis to "Submit" 				| X              |
| Add ViDGER download link 						|                |
| Add scatter plots to "Submit" 			| X              |
| Add static images 									|                |
| Sumbit ViDGER to Bioconductor 			|                |
| Add FAQ 														|  					     |
| System info 												| X					     |
| Fix axis issues (histogram)					| X					     |
| Add DEG overview                    | In progress... |