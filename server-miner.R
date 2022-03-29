

observe({
  
  inFileT <- input$datafileT;inFileB <- input$datafileB;inFileG <- input$datafileG
  if(is.null(inFileT) & is.null(inFileG)) {
    print("Please input a file contains your intersted genes")
  }
  if(is.null(inFileB)) {
    print("Using all annotated genes as background genes...")
  }
  if(!is.null(inFileT)&!is.null(inFileB)) {
    print("updating analyzed results")
    target <- inputDataReactive()$target
    back <- inputDataReactive()$back
  }
  if(!is.null(input$MS.dataset)){
    category <- input$MS.dataset
  }
  if(!is.null(inFileG)){
    print("updating analyzed results")
    gL <- inputDataReactive()$gL
  }

}
)


inputDataReactive <- reactive({
  
  print("inputting data")
  # Check if example selected, or if not then ask to upload a file.
  validate(
    need((input$lm_file_type=="exampleM")|(!is.null(input$datafileT)|(!is.null(input$datafileG))), 
         message = "Please input a file contains your intersted genes")
  )
  inFileT <- input$datafileT;inFileB <- input$datafileB;inFileG <- input$datafileG
  if ((!is.null(inFileT))&(!is.null(inFileB))) {
      target.list <- read.table(inFileT$datapath,header=F,stringsAsFactors = F)$V1
      back.list <- read.table(inFileB$datapath,header=F,stringsAsFactors = F)$V1
      print('uploaded gene list')
      return(list('target'=target.list));return(list('back'=back.list))
    } else if(!is.null(inFileG)){
      gL <- read.table(inFileG$datapath,header=F,stringsAsFactors = F)
      return(list('gL'=gL))
    } else  {
      return(NULL)
    }
  if(!is.null(input$MS.dataset)){
    return(list('category'=input$MS.dataset))
  }
})

# check if a file has been uploaded and create output variable to report this
output$fileUploaded <- reactive({
  return(!is.null(inputDataReactive()))
})
outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

# after the data is uploaded or example data is selected, analyze the data
analyzeDataReactive <- 
  eventReactive(input$upload_miner,
                ignoreNULL = FALSE, {
                  withProgress(message = "Analyzing data, please wait",{
                    if(input$lm_file_type=="exampleM"){
                      
                    } else {
                      
                    }
                  }
                  )
                }
  )

observeEvent(input$upload_miner, ({
  updateCollapse(session,id =  "input_miner_panel", open="miner_panel",
                 style = list("miner_panel" = "success",
                              "data_miner_panel"="primary"))
}))

output$mineroutput <- renderDataTable({
  print("output$mineroutput")
  getresults <- analyzeDataReactive()
  res = getresults$results
  datatable(res,rowname=F,width=10,class='compact hover',options = list(autoWidth  = T,
            initComplete = JS(
            "function(settings, json) {",
            "$('body').css({'font-family': 'Arial'});",
            "}"),
            initComplete = JS(
              "function(settings, json) {",
              "$(this.api().table().header()).css({'line-height': '50%';'width': '10px'});",
              "}"),
            columnDefs = list(list(targets= c(0:-1) ,render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 20 ?",
            "'<span title=\"' + data + '\">' +data.substr(0, 20) + '...</span>' : data;", "}"))),pageLength = 10, server = T)) %>% 
    formatStyle(c('p.adjust','pvalue','qvalues'),color = styleInterval(c(0.01, 0.05), c('#FF6B6B', '#FDD2BF', 'grey')))

})

# Download analyzed data

output$downloadResults_CSV <- downloadHandler(
  filename = paste0("EAGLE_results_",Sys.Date(),".csv"),
  content = function(file) {
    write.csv(analyzeDataReactive()$results, file, row.names=FALSE)})

output$downloadResults_RData <- downloadHandler(
  filename= paste0("EAGLE_results_",Sys.Date(),".RData"),
  content=function(file){
    start_list = analyzeDataReactive()
    save(start_list,file=file)
  })



