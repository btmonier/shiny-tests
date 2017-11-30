#------------------------------------------------------------------------------
# Title:  Shiny Test 11 - Crosstalk - User Interface
# Author: Brandon Monier (brandon.monier@sdstate.edu)
# Date:   11.24.17
#------------------------------------------------------------------------------

# Sources ----
source("tabs.R")


# User interface ----
vidgerUI <- fluidPage(
  	tabsetPanel(
  		tab.welcome,
  	  tab.submit,
  	  tab.deg,
  	  tab.heat,
  	  tab.help,
  	  tab.faq
  	)
)