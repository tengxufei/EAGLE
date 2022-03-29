

options(shiny.maxRequestSize = 100*1024^2)
#options(repos = BiocInstaller::biocinstallRepos()) # use setRepositories() 1 2

source("functions.R")
source("packages.R")
print(sessionInfo())

shinyServer(function(input, output, session) {
  ## Server functions are divided by tab
  ## 
  source("server-analysis.R",local = TRUE)
  #source("server-GSEA.R",local = TRUE)
  source("server-visulization.R",local = TRUE)
  #source("server-data.R",local = TRUE)
  source("server-geneanno.R",local = TRUE)
  #source("server-miner.R",local = TRUE)
  
})