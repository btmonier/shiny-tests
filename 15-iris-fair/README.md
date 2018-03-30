---
 Title:         IRIS Documentation
 Author:        Brandon Monier
 Created:       2018-03-30 at 09:35:53
 Last Modified: 2018-03-30 at 09:37:16
---

# IRIS - FAIR principle implementation

## About
This application is designed to test the functionality of outputting 
metadata to the GEO repository on NCBI.

## Get local app
Copy and paste this into your R environment:

``` r
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("shiny-tests", "btmonier", subdir = "15-iris-fair")
```