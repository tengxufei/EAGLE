tabPanel("Documentation",
         fluidRow(
           column(3,wellPanel(
             h4("Instructions"),
             h4("Enrichment Analysis",href = "#analysis"),
             a("Input data format", href = "#input"), br(),
             a("Functional data sets", href="#sets"), br(),
             a("Functional enrichment results", href="#result"), br(),
             h4("Visualization",href = "#visualization"),
             h4("Gene Annotation",href = "#anno"),
           )
           ),#column
           column(9,
                  includeMarkdown("instructions/Instructions.md"))
         ))