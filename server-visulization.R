


observe({
  print("server-analysisres-update")
  output$itemChoices <- renderUI({mydata = as.character(analyzeDataReactive()$plotResults$ID);selectInput('itemChoices', 'Select GeneSets to show', as.character(mydata),multiple = TRUE)})
  })


observe({
  
  print("drawing bar plot")
  
  data_analyzed = analyzeDataReactive()
  data_results = data_analyzed$results
  
  output$barplot <- renderPlotly({
    withProgress(message = "Drawing bar plot of selected items, please wait",
                 {
                   if (names(dev.cur()) != "null device") dev.off()
                   pdf(NULL)
                   p=barplot(res = analyzeDataReactive()$plotResults,
                                 absERcut = input$analysisres_enrich_ratio_cut,
                                 pvalcut = input$analysisres_pvalue_cut,
                                 fdrcut = input$analysisres_fdr_cut,
                                 color_low = input$barcolor_low,
                                 color_mid = input$barcolor_mid,
                                 color_high = input$barcolor_high,
                                 type=input$analysis_method)
                 }) 
    
  }) 
})

output$dotplot <- renderPlotly({
  withProgress(message = "Drawing bar plot of selected items, please wait",
               {
                 if (names(dev.cur()) != "null device") dev.off()
                 pdf(NULL)
                 p=dotplot(res = analyzeDataReactive()$plotResults,
                           absERcut = input$dot_enrich_ratio_cut,
                           pvalcut = input$dot_pvalue_cut,
                           fdrcut = input$dot_fdr_cut,
                           color_low = input$dotcolor_low,
                           color_mid = input$dotcolor_mid,
                           color_high = input$dotcolor_high,
                           type=input$analysis_method)
               }) 
  
}) 


observe({
  
  print("drawing GSEA plot")
  
  data_analyzed = analyzeDataReactive()
  GSEA_result <- analyzeDataReactive()$plotResults
  output$gseaWarning <- renderText({
    "Please select one or more GeneSets to show"
  })
  output$GSEAplot <- renderPlot({
    
    withProgress(message = "Drawing GSEA plot of selected items, please wait",
                 {
                   if(is.null(input$itemChoices)){
                     renderPrint("jj")
                   } else {
                     GSEA.plot(gseaResult = analyzeDataReactive()$plotResults,geneSetID = input$itemChoices,anno.cols = decode.col()[['choices']][[input$GSEAcolor]][1:length(input$itemChoices)])
                   }
                 }) 
    
  }) 
})

output$GSEAplot_PDF <- downloadHandler(
  filename=paste0("GSEAplot_",Sys.Date(),".pdf"),
  content = function(file) {
    ggsave(file=file,plot=GSEA.plot(gseaResult = analyzeDataReactive()$plotResults,geneSetID = input$itemChoices),width=11,height=7)
    })

output$networkplot_PDF <- downloadHandler(
  filename=paste0("networkplot_",Sys.Date(),".pdf"),
  content = function(file) {
    ggsave(file=file,plot=network.plot(res = analyzeDataReactive()$plotResults,
                                       absERcut = input$n_enrich_ratio_cut,
                                       pvalcut = input$n_pvalue_cut,
                                       fdrcut = input$n_fdr_cut,
                                       gene.list=analyzeDataReactive()$gene.list,
                                       color_low = input$ncolor_low,
                                       color_mid = input$ncolor_mid,
                                       color_high = input$ncolor_high,
                                       type=input$analysis_method),width=11,height=7)
  })

observe({
  
  print("drawing network")
  
  output$network_plot <- renderPlot({ 
    withProgress(message = "Drawing network, please wait",{
      #if (names(dev.cur()) != "null device") dev.off()
      #pdf(NULL)
      
      network.plot(res = analyzeDataReactive()$plotResults,
                absERcut = input$n_enrich_ratio_cut,
                pvalcut = input$n_pvalue_cut,
                fdrcut = input$n_fdr_cut,
                gene.list=analyzeDataReactive()$gene.list,
                color_low = input$ncolor_low,
                color_mid = input$ncolor_mid,
                color_high = input$ncolor_high,
                type=input$analysis_method)
      
    })
  })
  
  
  
  #}
})





