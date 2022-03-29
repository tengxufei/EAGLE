
source("functions.R")
source("packages.R")
source("subpathway.R")
library(shiny)

customHeaderPanel <- function(title,windowTitle=title){
  tagList(
    tags$head(
      tags$title(windowTitle),
      tags$link(rel="stylesheet", type="text/css",
                href="eagle.css"),
  
      tags$h1(a(href="eagle.css"))
    )
  )
}

tagList(
  tags$head(
    tags$link(HTML('<link rel="icon" type="image/png" sizes="144x144" href="eagle.png">')),
    tags$style(HTML(" .shiny-output-error-validation {color: darkred; } ")),
    tags$style(".mybuttonclass{background-color:#CD0000;} .mybuttonclass{color: #fff;} .mybuttonclass{border-color: #9E0000;}")
  ),
  navbarPage(
    
    theme = "eagle.css",
    title = "EAGLE: Enrichment analysis of genes in functional dataset",
    source("ui-doc.R",local=TRUE)$value,
    source("ui-analysis.R",local=TRUE)$value,
    source("ui-visulization.R",local=TRUE)$value,
    source("ui-geneanno.R",local=TRUE)$value,
    #source("ui-miner.R",local=TRUE)$value,
    ## ============================================================================ ##
    ## DOWNLOAD TAB
    ## ============================================================================ ##   
    #source("ui-tab-help.R",local=TRUE)$value,
    footer=p(hr(),p("Maintained by ***(Hidden due to blind review of thesis)", align="center",width=4),
             p(("Copyright (C) 2022, code licensed under GPLv3"),align="center",width=4),
             p(("Github:"),a("https://github.com/tengxufei/EAGLE",href="https://github.com/tengxufei/EAGLE"),align="center",width=4)
    ),
    
    tags$head(includeScript("a.js"))
  ) #end navbarpage
) #end taglist

