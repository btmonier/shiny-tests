# shiny-tests

## Overview

This is a collection of Shiny tests to determine the interactivity of the `vidger` package.

To run the examples locally, you can install the shiny package in R, and use the function `runGithub()`. For example, to run `13-vidger-final`:

``` r
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("shiny-tests", "btmonier", subdir = "14-iris-me")
```


## Fix list

```
Last updated: 2018-02-02 10:06:31 CST
```

| Task                                | Completed?     |
|-------------------------------------|----------------|
| Add clustering, heatmap to "Submit" | X              |
| Add DEG Analysis to "Submit" 		    | X              |
| Add scatter plots to "Submit"       | X              |
| Add static images 				          | X              |
| Sumbit ViDGER to Bioconductor       | In progress... |
| Add FAQ with ViDGER down link       | X              |
| System info                         | X              |
| Fix axis issues (histogram)         | X              |
| Add DEG overview                    | X              |
| Get link available                  | X              |
| Submit and QC split                 | X              |
| About us tab                        | X              |
