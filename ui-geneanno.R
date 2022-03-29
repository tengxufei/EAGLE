

tabPanel("Gene Annotation", 
         fluidRow(column(3,wellPanel(
           radioButtons('ga_file_type','Use example gene list or upload your own list',
                        choices = c('Upload Gene List'="uploadGL",'Example Genes'="exampleGL"),selected = "exampleGL"),
           radioButtons('annotation_method','Choose one method for gene annotation',
                        choices=c('Functional Gene Sets'="pathway",'Literature Mining'="literature"),selected = "pathway"),
           conditionalPanel(condition="input.ga_file_type=='uploadGL'",
                            fileInput('datafileGe', 'Choose file containing target gene list (.txt)',accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),),
           actionButton("upload_genelist","Submit Data",
                        style="color: #fff; background-color: #CD0000; border-color: #9E0000")
         )
         ),
         column(9,
                bsCollapse(id="geneanno_panel",open="dataGA_panel",multiple = FALSE,
                           bsCollapsePanel(title="Gene Annotation Results:",
                                           value="ga_panel",
                                           downloadButton('downloadResults_CSV2','Save Results as CSV File'),
                                           downloadButton('downloadResults_RData2',
                                                          'Save Results as .RData File for Future Analysis',
                                                          class="mybuttonclass"),
                                           DT::dataTableOutput('gaoutput'),
                                           tags$head(tags$style(".mybuttonclass{background-color:#CD0000;} .mybuttonclass{color: #fff;} .mybuttonclass{border-color: #9E0000;}"))
                           )
                )#bscollapse
         )#column
         )#fluidrow
)#tabpanel
