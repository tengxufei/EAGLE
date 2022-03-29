

observe({
  
  inFileGe <- input$datafileGe;
  if(is.null(inFileGe)) {
    print("Please input a file contains your intersted genes")
  } else {
    print("updating analyzed results")
  }
}
)

inputDataReactive.geneanno <- reactive({
  
  print("inputting data")
  # Check if example selected, or if not then ask to upload a file.
  validate(
    need((input$ga_file_type=="exampleGL")|(!is.null(input$datafileGe)), 
         message = "Please input a file contains your intersted genes")
  )
  inFileGe <- input$datafileGe
  if ((!is.null(inFileGe))) {
      genes <- read.table(inFileGe$datapath,header=F,stringsAsFactors = F)$V1
      print('uploaded gene list')
      return(list('genes'=genes));
    } else  {
      return(NULL)
    }
})

# check if a file has been uploaded and create output variable to report this
output$fileUploaded2 <- reactive({
  return(!is.null(inputDataReactive.geneanno()))
})
outputOptions(output, 'fileUploaded2', suspendWhenHidden=FALSE)

# after the data is uploaded or example data is selected, analyze the data
analyzeDataReactive.geneanno <- 
  eventReactive(input$upload_genelist,
                ignoreNULL = FALSE, {
                  withProgress(message = "Analyzing data, please wait",{
                    if(input$ga_file_type=="exampleGL"){
                        genes <- load_existing_rdata('data/exampleGenes.RData')
                    } else {
                        genes <- read.table(input$datafileGe$datapath,header=F,stringsAsFactors = F)$V1
                    }
                    res <- gene.anno(gene.list=genes,annotype=input$annotation_method)
                    return(list("results"=res,"plotResults"=res))
                    }
                  )
                }
  )

observeEvent(input$upload_genelist, ({
  updateCollapse(session,id =  "geneanno_panel", open="ga_panel",
                 style = list("ga_panel" = "success",
                              "dataGA_panel"="primary"))
}))

output$gaoutput <- DT::renderDataTable({
  print("output$gaoutput")
  getresultsGA <- analyzeDataReactive.geneanno()
  resGA = getresultsGA$results
  DT::datatable(resGA,
                escape = FALSE,rowname=F,width=10,class='compact hover',options = list(autoWidth  = T,
            initComplete = JS(
            "function(settings, json) {",
            "$('body').css({'font-family': 'Arial'});",
            "}"),
            initComplete = JS(
              "function(settings, json) {",
              "$(this.api().table().header()).css({'line-height': '50%';'width': '10px'});",
              "}"),
            columnDefs = list(list(targets= c(1:5,-2,-1) ,render = JS(
            "function(data, type, row, meta) {",
            "return type === 'display' && data.length > 20 ?",
            "'<span title=\"' + data + '\">' +data.substr(0, 20) + '...</span>' : data;", "}"))),pageLength = 10, server = T)) 
})

# Download analyzed data

output$downloadResults_CSV2 <- downloadHandler(
  filename = paste0("EAGLE_results_",Sys.Date(),".csv"),
  content = function(file) {
    write.csv(analyzeDataReactive.geneanno()$results, file, row.names=FALSE)})

output$downloadResults_RData2 <- downloadHandler(
  filename= paste0("EAGLE_results_",Sys.Date(),".RData"),
  content=function(file){
    start_list = analyzeDataReactive.geneanno()
    save(start_list,file=file)
  })



