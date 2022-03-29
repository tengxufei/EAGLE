

tabPanel("Literature knowledge", 
         fluidRow(column(3,wellPanel(
           radioButtons('lm_file_type','Use example file or upload your own file',
                        choices = c('Upload PubMed Abstract File'="uploadM",'Example Data'="exampleM"),selected = "exampleM"),
           conditionalPanel(condition="input.lm_file_type=='uploadM'",
                            
                            
           ),
           conditionalPanel(condition="input.lm_file_type=='exampleM'",
                            textInput(inputId = "bp",
                                         "Biological phrase", value = 25)),
           

           actionButton("upload_miner","Submit Data",
                        style="color: #fff; background-color: #CD0000; border-color: #9E0000")
         )
         ),
         ## ==================================================================================== ##
         ## Right hand column shows data input DT and data analysis result DT
         ## ==================================================================================== ##
         column(9,
                bsCollapse(id="input_miner_panel",open="data_miner_panel",multiple = FALSE,
                           bsCollapsePanel(title="Enrichment Analysis Results:",
                                           value="miner_panel",
                                           downloadButton('downloadResults_CSV','Save Results as CSV File'),
                                           downloadButton('downloadResults_RData',
                                                          'Save Results as .RData File for Future Analysis',
                                                          class="mybuttonclass"),
                                           dataTableOutput('mineroutput'),
                                           tags$head(tags$style(".mybuttonclass{background-color:#CD0000;} .mybuttonclass{color: #fff;} .mybuttonclass{border-color: #9E0000;}"))
                           )
                )#bscollapse
         )#column
         )#fluidrow
)#tabpanel
