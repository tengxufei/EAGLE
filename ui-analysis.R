tabPanel("Enrichment Analysis", 
         fluidRow(column(3,wellPanel(
           radioButtons('data_file_type','Use example file or upload your own file',
                        choices = c('Upload Data'="upload",'Example Data'="exampleGenes"),selected = "exampleGenes"),
           radioButtons('analysis_method','Choose one enrichment analysis method',
                        choices=c('Over-Representation Analysis'="hyperG",'GSEA analysis'="gsea",'Subpathway analysis'="subpath"),selected = "gsea"),
           conditionalPanel(condition="input.data_file_type=='upload'",
                            conditionalPanel(condition="input.analysis_method!='gsea'",
                                             fileInput('datafileT', 'Choose file containing target gene list (.txt)',
                                                       accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                                             fileInput('datafileB', 'Choose file containing background gene list (.txt)',
                                                       accept=c('text/csv', 'text/comma-separated-values,text/plain','.txt')),               
                            ),
                            conditionalPanel(condition="input.analysis_method=='gsea'",
                                             fileInput('datafileG', 'Choose file containing target gene matrix (.txt)',
                                                       accept=c('text/csv','text/comma-separated-values,text/plain','.txt')),
                            ),
                            
           ),
           conditionalPanel(condition="input.data_file_type=='exampleGenes'",),
           conditionalPanel(condition="input.analysis_method!='subpath'",
             selectInput("database", label = h4("Collected Gene Sets"), 
                         choices = list("m6A regulator"="MR","ncRNA"="LA","MSigDB"="MS","Kinase"="KI")),
             conditionalPanel(condition="input.database=='MS'",
                              selectInput("MS.dataset", label = h4("MSigDB dataset"), 
                                          choices = list("H (hallmark gene sets)"="H","C1 (positional gene sets)"="C1","C2 (curated gene sets)"="C2","C3 (regulatory target gene sets)"="C3","C4 (computational gene sets)"="C4","C5 (ontology gene sets)"="C5","C6 (oncogenic signature gene sets)"="C6","C7 (immunologic signature gene sets)"="C7","C8 (cell type signature gene sets)"="C8")),
             ),),
           conditionalPanel(condition="input.analysis_method=='subpath'",
                            selectInput("database2", label = h4("Collected Gene Sets"),
                              choices = list("Kyoto Encyclopedia of Genes and Genomes"="KEGG")),
                            sliderInput(inputId="k.num", "Number of kCliques", 1, 10, 4)),

           actionButton("upload_data","Submit Data",
                        style="color: #fff; background-color: #CD0000; border-color: #9E0000")
         )
         ),
         ## ==================================================================================== ##
         ## Right hand column shows data input DT and data analysis result DT
         ## ==================================================================================== ##
         column(9,
                bsCollapse(id="input_collapse_panel",open="data_panel",multiple = FALSE,
                           bsCollapsePanel(title="Enrichment Analysis Results:",
                                           value="analysis_panel",
                                           downloadButton('downloadResults_CSV','Save Results as CSV File'),
                                           downloadButton('downloadResults_RData',
                                                          'Save Results as .RData File for Future Analysis',
                                                          class="mybuttonclass"),
                                           DT::dataTableOutput('analysisoutput'),
                                           tags$head(tags$style(".mybuttonclass{background-color:#CD0000;} .mybuttonclass{color: #fff;} .mybuttonclass{border-color: #9E0000;}"))
                           )
                )#bscollapse
         )#column
         )#fluidrow
)#tabpanel
