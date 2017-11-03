#------------------------------------------------------------------------------
# Title:  Shiny Test 07 - Modularity - User Interface
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   10.26.17
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