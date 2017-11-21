#------------------------------------------------------------------------------
# Title:  Shiny Test 08 - Crosstalk - User Interface
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.03.17
#------------------------------------------------------------------------------

# Sources ----
source("tabs.R")


# User interface ----
vidgerUI <- fluidPage(
  	tabsetPanel(
  		tab.welcome,
  	  	tab.submit,
  	  	tab.deg,
  	  	tab.help
  	)
)