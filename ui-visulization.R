## tengxufei@gmail.com

tabPanel("Enrichment Visulization",  
         fluidRow(column(3,wellPanel(
           conditionalPanel("input.analysisres_tabset=='Bar plot'",
                            #selectizeInput("analysisres_type", label="Select Enriched Results for Bar Plot",
                                           #choices=c('Hypergeometric test'="hyperG",'GSEA analysis'="gsea"),selected = "analysis_method"),
                            numericInput("analysisres_enrich_ratio_cut",
                                         label="Choose Enrich Ratio Threshold",min= 0, max=1,value=0.1),
                            numericInput("analysisres_pvalue_cut",
                                         label="Choose P-value Threshold",min=0,max=1,value=0.05),
                            numericInput("analysisres_fdr_cut",
                                         label="Choose adjusted P-value Threshold",min=0,max=1,value=0.05),
                            colourInput("barcolor_low",label="Select Color - Low Values","#51adcf",
                                        showColour = "background"),
                            colourInput("barcolor_mid",label="Select Color - Middle Values","grey",
                                        showColour = "background"),
                            colourInput("barcolor_high",label="Select Color - High Values","#d3de32",
                                        showColour = "background")
           ),
           conditionalPanel("input.analysisres_tabset=='Dot plot'",
                            #selectizeInput("analysisres_type", label="Select Enriched Results for Bar Plot",
                            #choices=c('Hypergeometric test'="hyperG",'GSEA analysis'="gsea"),selected = "analysis_method"),
                            numericInput("dot_enrich_ratio_cut",
                                         label="Choose Enrich Ratio Threshold",min= 0, max=1,value=0.1),
                            numericInput("dot_pvalue_cut",
                                         label="Choose P-value Threshold",min=0,max=1,value=0.05),
                            numericInput("dot_fdr_cut",
                                         label="Choose adjusted P-value Threshold",min=0,max=1,value=0.05),
                            colourInput("dotcolor_low",label="Select Color - Low Values","#51adcf",
                                        showColour = "background"),
                            colourInput("dotcolor_mid",label="Select Color - Middle Values","grey",
                                        showColour = "background"),
                            colourInput("dotcolor_high",label="Select Color - High Values","#d3de32",
                                        showColour = "background")
           ),
           conditionalPanel("input.analysisres_tabset=='GSEA plot'",
                            uiOutput('itemChoices'),
                            palettePicker(
                              inputId = "GSEAcolor", label = "Select palette for GSEA items", 
                              choices = list(
                                  "viridis" = viridis_pal(option = "viridis")(8),
                                  "magma" = viridis_pal(option = "magma")(8),
                                  "inferno" = viridis_pal(option = "inferno")(8),
                                  "plasma" = viridis_pal(option = "plasma")(8),
                                  "cividis" = viridis_pal(option = "cividis")(8),
                                  "Blues" = brewer_pal(palette = "Blues")(8),
                                  "Reds" = brewer_pal(palette = "Reds")(8),
                                  "Paired" = brewer_pal(palette = "Paired")(8),
                                  "Set1" = brewer_pal(palette = "Set1")(8),
                                  "Set2" = brewer_pal(palette = "Set2")(8),
                                  "Set3" = brewer_pal(palette = "Set3")(8)
                              ), 
                              textColor = c(
                                rep("white", 5), rep("black", 4) 
                              )
                            ),
                            radioButtons("pvalue_table",label="whether add a pvalue table",
                                         choices=c("Yes","No"),selected = "No")
                            
           ),
           conditionalPanel("input.analysisres_tabset=='DIY plot by esquisse'",
                            numericInput(inputId = "n",
                                         "Sample size", value = 25)
           ),
           conditionalPanel("input.analysisres_tabset=='Network plot'",
                            #selectizeInput("n_type", label="Select Enriched Results for cnet Plot",
                             #               choices=c('Hypergeometric test'="hyperG",'GSEA analysis'="gsea"),selected = "analysis_method"),
                            numericInput("n_enrich_ratio_cut",
                                         label="Choose Enrich Ratio Threshold",min= 0, max=1,value=0),
                            numericInput("n_pvalue_cut",
                                         label="Choose P-value Threshold",min=0,max=1,value=1),
                            numericInput("n_fdr_cut",
                                         label="Choose adjusted P-value Threshold",min=0,max=1,value=1),
                            colourInput("ncolor_low",label="Select Color - Low Values","#51adcf",
                                        showColour = "background"),
                            colourInput("ncolor_mid",label="Select Color - Middle Values","grey",
                                        showColour = "background"),
                            colourInput("ncolor_high",label="Select Color - High Values","#d3de32",
                                        showColour = "background")
                            
           )
           
         )#,
         ),#column
         column(9,
                tabsetPanel(id="analysisres_tabset",
                            tabPanel(title="Bar plot",
                                     plotlyOutput("barplot",height=600)
                            ),
                            tabPanel(title="Dot plot",
                                     plotlyOutput("dotplot",height=600)
                            ),
                            tabPanel(title="Network plot",
                                     downloadButton('networkplot_PDF','Save Plot as PDF File'),
                                     plotOutput(outputId = "network_plot")
                            ),
                            tabPanel(title="GSEA plot",
                                     downloadButton('GSEAplot_PDF','Save Plot as PDF File'),
                                     if(is.null("input.itemChoices")) {
                                       textOutput('gseaWarning',container = 'h1')
                                     } else {
                                       plotOutput("GSEAplot",height=600)
                                     }
                                     
                            )
                            
                          
                )#tabsetPanel
         )#column
         )#fluidrow
) #END tabPanel