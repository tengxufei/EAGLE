getAnn<-function(geneList,background=getDefaultBackground(),
   order="pvalue",decreasing=FALSE,graphList=getDefaultGraph()){
      if(typeof(geneList)!="character"){
          print("warning: your geneList must be 'character' vector. Because the type of your current geneList is not correct, it has been conveted arbitrarily using the function as.character().")
          geneList<-as.character(geneList)
          }
      if(!exists("ke2g")) initialize_ke2g()
      graphList<-graphList[sapply(graphList,function(x) length(x)>0)]     
      keggpathid2name<-get("keggpathid2name",envir=ke2g)

      annList<-list()
      for(i in 1:length(graphList)){
            ann<-list(pathwayName="not known",annGeneList=character(),annBgGeneList=character(),annGeneNumber=0,
                      annBgNumber=0,geneNumber=0,bgNumber=0,pvalue=1,qvalue=1)
            if(class(graphList[[i]])=="character"){
                  graphGeneList<-getGeneFromPathway(graphList[[i]])
                                  pathwayName<-getPNameFromPId(graphList[[i]])
            }
            else if(class(graphList[[i]])=="graphNEL"){
                  graphGeneList<-gsub("hsa:","",getKGeneFromKO(nodes(graphList[[i]])))
                  pathwayName<-keggpathid2name[[substring(names(graphList)[i],6,10)]]                             
            }            
            annotatedGeneList<-intersect(graphGeneList,geneList)
            annotatedBackgroundList<-intersect(graphGeneList,background)

            #pathwayName<-keggpathid2name[names(graphList)[i]]
            if(length(pathwayName)!=0)
                ann$pathwayName<-pathwayName
            ann$annGeneList<-annotatedGeneList   
            ann$annBgGeneList<-annotatedBackgroundList
            ann$annGeneNumber<-length(annotatedGeneList)
            ann$annBgNumber<-length(annotatedBackgroundList)

            ann$geneNumber<-length(geneList)
            ann$bgNumber<-length(background)

            ann$pvalue<-1-phyper(ann$annGeneNumber-1,ann$annBgNumber,
                 ann$bgNumber-ann$annBgNumber,ann$geneNumber)
            
            annList[[i]]<-ann
      } 
qvalueList<-fdrtool(sapply(annList,function(x) x$pvalue), 
                   statistic="pvalue",plot=FALSE,verbose=FALSE)$qval
      for(i in 1:length(annList)){
            annList[[i]]$qvalue<-qvalueList[i]
      }
      names(annList)<-names(graphList)
      annList<-annList[sapply(annList,function(x) x$annGeneNumber>0)]
      annList<-annList[order(sapply(annList,function(x) x[[order]]),decreasing=decreasing)]

      return(annList)
}
getDefaultUndirectedGraph<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      uGraph<-get("uGraph",envir=ke2g)
      return(uGraph)       
}
getGeneFromPathway<-function(pathwayList){
      return(getGeneFromKGene(getKGeneFromPathway(pathwayList)))
}
getKGeneFromPathway<-function(pathwayList){
	  pathwayList<-as.character(pathwayList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2path<-get("gene2path",envir=ke2g)  
      keggGeneList<-unique(as.character(gene2path[as.character(gene2path[,2]) %in% paste("path:",getOrgAndIdType()[1],substring(pathwayList,6),sep=""),1]))
      return(keggGeneList)
}
getOrgAndIdType<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      orgAndIdType<-get("orgAndIdType",envir=ke2g)
      return(orgAndIdType)    
}
getPNameFromPId<-function(PIdList){
	  PIdList<-as.character(PIdList)
      if(!exists("ke2g")) initialize_ke2g()
	  map_title<-get("map_title",envir=ke2g)
      PNameList<-unique(as.character(map_title[as.character(map_title[,1]) %in% substring(PIdList,6),2]))
      return(PNameList)
}

getKcSubGraph<-function(k=4,graphList=getDefaultKOUndirectedGraph()){
      if(k<1)
            stop("k can't <1 ")
      graphList<-graphList[sapply(graphList,function(x) length(x)>0)]
      subGraphIndex<-0
      subGraph<-list()
      subNames<-character()
      for(i in 1:length(graphList)){
         kc<-kCliques(graphList[[i]])
         if(length(kc)>0){
            if(k<=length(kc)){
                  kkc<-kc[[k]]
            }
            else{
                  kkc<-kc[[length(kc)]]
            }

            for(j in 1:length(kkc)){
                  subGraphIndex<-subGraphIndex+1
                  subGraph[subGraphIndex]<-subGraph(kkc[[j]],graphList[[i]])
                  subNames[subGraphIndex]<-paste(names(graphList)[i],j,sep="_")
            }
          }
      }
      names(subGraph)<-subNames
      return(subGraph)
}

getDefaultKOUndirectedGraph<-function(){
      if(!exists("ke2g")) initialize_ke2g()
      KOuGraph<-get("KOuGraph",envir=ke2g)
      return(KOuGraph)       
}

getDefaultBackground<-function(){
      if(!exists("ke2g")) initialize_ke2g()
	  keggGene2gene<-get("keggGene2gene",envir=ke2g) 
	  background<-unique(as.character(keggGene2gene[,2]))
	  newBackground<-sapply(strsplit(background,":"),function(x) x[2])
      return(newBackground)
}

getGeneFromEnzyme<-function(enzymeList){
      return(getGeneFromKGene(getKGeneFromEnzyme(enzymeList)))
}


getKGeneFromEnzyme<-function(enzymeList){
    enzymeList<-as.character(enzymeList)
    if(!exists("ke2g")) initialize_ke2g()
    gene2ec<-get("gene2ec",envir=ke2g)
    keggGeneList<-unique(as.character(gene2ec[as.character(gene2ec[,2]) %in% enzymeList,1]))
    return(keggGeneList)
}

getDefaultUndirectedGraph<-function(){
      uGraph<-get("uGraph",envir=ke2g)
      return(uGraph)       
}

getAexample<-function(k=1000){
   if(!exists("ke2g")) initialize_ke2g()
   gene2path<-get("gene2path",envir=ke2g)
   allGene1<-getGeneFromKGene(as.character(gene2path[,1]))   
   if(k<=length(allGene1))
   geneList<-allGene1[1:k]
   else
   geneList<-allGene1

   return(geneList)
}

getGeneFromKGene<-function(keggGeneList){
	  keggGeneList<-as.character(keggGeneList)
      if(!exists("ke2g")) initialize_ke2g()
      keggGene2gene<-get("keggGene2gene",envir=ke2g)
      geneList<-unique(as.character(sapply(strsplit(as.character(keggGene2gene[as.character(keggGene2gene[,1]) %in% keggGeneList,2]),":"),function(x) return (x[2]))))
      return(geneList)
}
getKGeneFromKO<-function(KOList){
	  KOList<-as.character(KOList)
      if(!exists("ke2g")) initialize_ke2g()
	  gene2ko<-get("gene2ko",envir=ke2g)
      keggGeneList<-unique(as.character(gene2ko[as.character(gene2ko[,2]) %in% KOList,1]))
      return(keggGeneList)
}
