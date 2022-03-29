

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
    category <- inputDataReactive()$category
    ## it doesn't work, the "category <- input$MS.dataset" works
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
    need((input$data_file_type=="exampleGenes")|(!is.null(input$datafileT)|(!is.null(input$datafileG))), 
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
#  if(!is.null(input$MS.dataset)){
#    return(list('category'=input$MS.dataset))
#  }
})

# check if a file has been uploaded and create output variable to report this
output$fileUploaded <- reactive({
  return(!is.null(inputDataReactive()))
})
outputOptions(output, 'fileUploaded', suspendWhenHidden=FALSE)

# after the data is uploaded or example data is selected, analyze the data
analyzeDataReactive <- 
  eventReactive(input$upload_data,
                ignoreNULL = FALSE, {
                  withProgress(message = "Analyzing data, please wait",{
                    if(input$data_file_type=="exampleGenes"){
                      if(input$analysis_method=="hyperG") {
                        target <- load_existing_rdata('data/exampleGenes.RData')
                        back <- load_existing_rdata('data/exampleBack.RData')
                        category <- input$MS.dataset
                      } else if(input$analysis_method=="subpath") {
                        target <- load_existing_rdata('data/exampleGenes.RData')
                        back <- load_existing_rdata('data/exampleBack.RData')
                        k.num <- input$k.num
                      } else if(input$analysis_method=="gsea"){
                        gL <- load_existing_rdata('data/gL.RData')
                        if(!is.null(input$MS.dataset)){
                          category <- input$MS.dataset
                        }
                      }
                    } else {
                      if(input$analysis_method=="hyperG"){
                        target.list <- read.table(input$datafileT$datapath,header=F,stringsAsFactors = F)$V1
                        back.list <- read.table(input$datafileB$datapath,header=F,stringsAsFactors = F)$V1
                        category <- input$MS.dataset
                        } else if(input$analysis_method=="gsea"){
                        gL <- inputDataReactive()$gL
                        if(!is.null(input$MS.dataset)){
                          category <- input$MS.dataset
                        } 
                      } else if(input$analysis_method=="subpath"){
                        target.list <- read.table(input$datafileT$datapath,header=F,stringsAsFactors = F)$V1
                        back.list <- read.table(input$datafileB$datapath,header=F,stringsAsFactors = F)$V1
                        k.num <- input$k.num
                      }
                    }
                    if(input$analysis_method=="hyperG") {
                      if(input$database=="MR") {
                        myfiles <- paste("data/input_dataset/MR",list.files(path = "data/input_dataset/MR"),sep="/")
                        data_all <- lapply(myfiles, function(x){read.table(x,header = F,sep="\t",stringsAsFactors = F)})
                        names(data_all) <- gsub("data/input_dataset/MR/","",myfiles)
                        item <- do.call(rbind,lapply(c(1:length(data_all)),function(x){data.frame(V1=names(data_all)[x],V2=data_all[[x]][,1])}))
                        res <- fdr.based.hyperG.test(fun.df=item,gene.list=target,back.list=back)
                        return(list("results"=res,"plotResults"=res))
                      }else if(input$database=="MS") {
                        genesets = msigdbr(species = "Homo sapiens", category = category)
                        t2g = genesets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
                        res <- fdr.based.hyperG.test(fun.df=t2g,gene.list=target,back.list=back)
                        return(list("results"=res,"plotResults"=res))
                      }else if(input$database=="LA") {
                        nc <- read.table("data/input_dataset/ncRNA.inter.txt",header = T,stringsAsFactors = F,sep="\t")
                        res <- fdr.based.hyperG.test(fun.df=nc,gene.list=target,back.list=back)
                       # table(item$V1)<10
                        #增加min size 和max size的筛选
                        #还有一个问题是那些"excel基因"要弄好 
                        return(list("results"=res,"plotResults"=res))
                      }else if(input$database=="KI") {
                        ks <- read.table("data/input_dataset/kinase.substrate.txt",header = T,stringsAsFactors = F,sep="\t")
                        res <- fdr.based.hyperG.test(fun.df=ks,gene.list=target,back.list=back)
                        return(list("results"=res,"plotResults"=res))
                      }
                    } 
                    if(input$analysis_method=="gsea") {
                      if(input$database=="MR") {
                        myfiles <- paste("data/input_dataset/MR",list.files(path = "data/input_dataset/MR"),sep="/")
                        data_all <- lapply(myfiles, function(x){read.table(x,header = F,sep="\t",stringsAsFactors = F)})
                        names(data_all) <- gsub("data/input_dataset/MR/","",myfiles)
                        item <- do.call(rbind,lapply(c(1:length(data_all)),function(x){data.frame(V1=names(data_all)[x],V2=data_all[[x]][,1])}))
                        resG <- GSEA.analysis(category=item,gene.list=gL)
                        return(list("results"=resG[[1]],"plotResults"=resG[[2]],"gene.list"=gL))
                      }else if(input$database=="MS") {
                        resG <- GSEA.MSigDB(category=category,gene.list=gL)
                        return(list("results"=resG[[1]],"plotResults"=resG[[2]],"gene.list"=gL))
                      }else if(input$database=="LA") {
                        nc <- read.table("data/input_dataset/ncRNA.inter.txt",header = T,stringsAsFactors = F,sep="\t")
                        resG <- GSEA.analysis(category=nc,gene.list=gL)
                        return(list("results"=resG[[1]],"plotResults"=resG[[2]],"gene.list"=gL))
                      }else if(input$database=="KI") {
                        ks <- read.table("data/input_dataset/kinase.substrate.txt",header = T,stringsAsFactors = F,sep="\t")
                        resG <- GSEA.analysis(category=ks,gene.list=gL)
                        return(list("results"=resG[[1]],"plotResults"=resG[[2]],"gene.list"=gL))
                      }
                    }
                    if(input$analysis_method=="subpath") {
                      if(input$database2=="KEGG") {
                        resG <- subpath.KEGG(target.list=target,back.list=back,kNum=k.num)
                        #resG <- data.frame(data=1:10,value=1:10)
                        return(list("results"=resG,"plotResults"=resG))
                      }
                    }
                  }
                  )
                }
  )

observeEvent(input$upload_data, ({
  updateCollapse(session,id =  "input_collapse_panel", open="analysis_panel",
                 style = list("analysis_panel" = "success",
                              "data_panel"="primary"))
}))

output$analysisoutput <- DT::renderDataTable({
    print("output$analysisoutput")
    getresults <- analyzeDataReactive()
    res = getresults$results
    DT::datatable(res,rowname=F,class='compact hover',options = list(autoWidth  = T,scrollX=TRUE, scrollCollapse=TRUE,serverSide=FALSE,
                                                                          initComplete = JS(
                                                                            "function(settings, json) {",
                                                                            "$('body').css({'font-family': 'Arial'});",
                                                                            "}"),
                                                                          initComplete = JS(
                                                                            "function(settings, json) {",
                                                                            "$(this.api().table().header()).css({'line-height': '50%';'width': '10px'});",
                                                                            "}"),
                                                                          columnDefs = list(list(targets= c(0:5,-1) ,render = JS(
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



