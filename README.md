# shiny-tests

## Overview

This is a collection of Shiny tests to determine the interactivity of the `vidger` package.

To run the examples locally, you can install the shiny package in R, and use the function `runGithub()`. For example, to run `14-iris-me`:

``` r
if (!require('shiny')) install.packages("shiny")
shiny::runGitHub("shiny-tests", "btmonier", subdir = "14-iris-me")
```
