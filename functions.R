#### functions ~

load_existing_rdata <- function(rdata_filepath) {
  start_data <- load(rdata_filepath)
  example.result <-  get(start_data)
  return(example.result)
}

fdr.based.hyperG.test <- function(fun.df,gene.list,back.list) {
  
  
  fun.df <- fun.df[fun.df[,2] %in% back.list,]
  
  gene2item=tapply(fun.df[,1],as.factor(fun.df[,2]),function(x) unique(x))
  item2gene=tapply(fun.df[,2],as.factor(fun.df[,1]),function(x) unique(x))

  geneWithItem <- intersect(gene.list,names(gene2item)) ## annotated target genes
  n=length(geneWithItem);N=length(gene2item) ## n:target genes; N: back genes
  hyperG.test <- function(item,ix=item2gene,gwi=geneWithItem) {
    M=length(ix[[item]])  ## number of genes in this item
    exp.ratio=n*M/N  ## expected ratio of target genes in each item
    geneid=intersect(gwi,ix[[item]])
    k=length(geneid)
    odds.ratio=k/exp.ratio
    p=as.numeric(sprintf("%7.2e",phyper(k-1, M, N-M, n, lower.tail=F)))
    data.frame(Description=item,pvalue=p,odds.ratio=odds.ratio,enrich.ratio=as.numeric(sprintf("%7.2e",(k/M))),counts=k,geneID=paste(geneid,collapse = '/'))
  }
  res <- do.call(rbind,lapply(names(item2gene),function(x){hyperG.test(item=x)}))
  res$p.adjust <- as.numeric(sprintf("%7.2e",p.adjust(res$pvalue,"BH")))
  res$qvalues <- as.numeric(sprintf("%7.2e",qvalue(p=res$pvalue, lambda=0.05, pi0.method="bootstrap")$qvalue))
  res$odds.ratio <- as.numeric(sprintf("%7.2e",res$odds.ratio))
  #res$enrich.ratio <- sprintf("%7.2e",res$enrich.ratio)
  rownames(res) <- NULL
  res[,c(1,3:5,2,7:8,6)]
}

GSEA.MSigDB <- function(category,gene.list) {
  geneList <- gene.list[,2]
  names(geneList) <- gene.list[,1]
  geneList <- sort(geneList,decreasing = T)
  genesets = msigdbr(species = "Homo sapiens", category = category)
  t2g = genesets %>% dplyr::select(gs_name, gene_symbol) %>% as.data.frame()
  ret.mat <-  GSEA(geneList, TERM2GENE = t2g,pvalueCutoff = 1)
  rownames(ret.mat) <- NULL
  r2 <- ret.mat;ret.mat <- as.data.frame(ret.mat)
  ret.mat$enrichmentScore <- as.numeric(sprintf("%7.2e",ret.mat$enrichmentScore))
  ret.mat$NES <- as.numeric(sprintf("%7.2e",ret.mat$NES))
  ret.mat$pvalue <- as.numeric(sprintf("%7.2e",ret.mat$pvalue))
  ret.mat$p.adjust <- as.numeric(sprintf("%7.2e",ret.mat$p.adjust))
  ret.mat$qvalues <- as.numeric(sprintf("%7.2e",ret.mat$qvalues))
  list(ret.mat[,-1],r2)
}

GSEA.analysis <- function(category,gene.list) {
  geneList <- gene.list[,2]
  names(geneList) <- gene.list[,1]
  geneList <- sort(geneList,decreasing = T)
  ret.mat <- GSEA(geneList, TERM2GENE = category,pvalueCutoff = 1)
  r2 <- ret.mat;ret.mat <- as.data.frame(ret.mat)
  ret.mat$enrichmentScore <- as.numeric(sprintf("%7.2e",ret.mat$enrichmentScore))
  ret.mat$NES <- as.numeric(sprintf("%7.2e",ret.mat$NES))
  ret.mat$pvalue <- as.numeric(sprintf("%7.2e",ret.mat$pvalue))
  ret.mat$p.adjust <- as.numeric(sprintf("%7.2e",ret.mat$p.adjust))
  ret.mat$qvalues <- as.numeric(sprintf("%7.2e",ret.mat$qvalues))
  list(ret.mat[,-1],r2)
}

subpath.KEGG <- function(target.list,back.list,kNum){
  tmpG <- bitr(target.list, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
  bacG <- bitr(back.list, fromType = "SYMBOL", toType = "ENTREZID",OrgDb = "org.Hs.eg.db")$ENTREZID
  subG <- getKcSubGraph(k=kNum,graphList=getDefaultKOUndirectedGraph())
  ann2 <- getAnn(geneList=tmpG,background=bacG,order="pvalue",decreasing=FALSE,graphList=subG)
  ann2 <- as.data.frame(do.call(rbind,ann2))
  ann2$pathwayName <- unlist(ann2$pathwayName)
  ann2$marker <- unlist(ann2$marker)
  ann2$pvalue <- unlist(ann2$pvalue)
  ann2$p.adjust <- as.numeric(sprintf("%7.2e",p.adjust(ann2$pvalue,"BH")))
  ann2$annGeneList <- unlist(lapply(ann2$annGeneList,function(x){paste(x,collapse="/")}))
  ann2$annBgGeneList <- unlist(lapply(ann2$annBgGeneList,function(x){paste(x,collapse="/")}))
  test <- lapply(ann2$annGeneList,function(x){paste(bitr(strsplit(x,split = "/")[[1]], fromType = "ENTREZID", toType = "SYMBOL",OrgDb = "org.Hs.eg.db")$SYMBOL,collapse = "/")})
  ann2$annGeneList <- unlist(test)
  
  test <- lapply(ann2$annBgGeneList,function(x){paste(bitr(strsplit(x,split = "/")[[1]], fromType = "ENTREZID", toType = "SYMBOL",OrgDb = "org.Hs.eg.db")$SYMBOL,collapse = "/")})
  ann2$annBgGeneList <- unlist(test)
  ann2$annGeneNumber <- unlist(ann2$annGeneNumber)
  ann2$subpathwayID <- rownames(ann2)
  ann2$annBgNumber <- unlist(ann2$annBgNumber)
  ann2$geneNumber <- unlist(ann2$geneNumber)
  ann2$bgNumber <- unlist(ann2$bgNumber)
  colnames(ann2) <- gsub('qvalue','qvalues',colnames(ann2))
  ann2$qvalues <- unlist(ann2$qvalues)
  ann2
}

gene.anno <- function(gene.list,annotype){
  
  if(annotype=="pathway"){
    myfiles <- paste("data/input_dataset/MR",list.files(path = "data/input_dataset/MR"),sep="/")
    data_all <- lapply(myfiles, function(x){read.table(x,header = F,sep="\t",stringsAsFactors = F)})
    names(data_all) <- gsub("data/input_dataset/MR/","",myfiles)
    item <- do.call(rbind,lapply(c(1:length(data_all)),function(x){data.frame(V1=names(data_all)[x],V2=data_all[[x]][,1])}))
    
    kegg <- read.table("data/KEGG_hsa.txt",header = F,stringsAsFactors = F,sep="\t")
    msigdb <- as.data.frame(fread("data/msigdb_hsa.txt",header = F,stringsAsFactors = F,sep="\t"))
    
    load('data/n3.Rdata')
    n3 <- n3[n3$hgnc_symbol %in% gene.list,c(1,3,2)];colnames(n3) <- c("symbol","description","type")
    n3[is.na(n3)] <- ""
    n3$m6ARegulator <- unlist(lapply(n3$symbol,function(x){paste(item[item$V2==x,1],collapse = ";")}))
    n3$KEGGPathway <- gsub(" ;",";",unlist(lapply(n3$symbol,function(x){paste(kegg[kegg$V2==x,1],collapse = ";")})))
    n3$MSigDB <- tolower(unlist(lapply(n3$symbol,function(x){paste(msigdb[msigdb$V2==x,1],collapse = ";")})))
    #n3$hgnc_symbol = sapply(n3$hgnc_symbol,function(x){toString(tags$a(href=paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",x), x))})
    #path <- paste0("https://www.genecards.org/cgi-bin/carddisp.pl?gene=",n3$hgnc_symbol)
    #n3$hgnc_symbol <- paste0('<a href=', path, '>', n3$hgnc_symbol,'</a>' )
  }
  if(annotype=="literature"){
    load('data/lm.Rdata')
    n3 <- lm[lm$symbol %in% gene.list,]
  }
  n3
} 



barplot <- function(res, absERcut=0.05,pvalcut=0.05,fdrcut=0.05,color_low="#51adcf",
                        color_mid="grey",color_high="#d3de32",type) {
  
  if(type=="gsea") {
    res <- res[res$pvalue<pvalcut&res$p.adjust<fdrcut,]
    res <- res[order(res$enrichmentScore),]
    res$ID <- factor(res$Description,levels=res$Description,ordered = T)
    p <- ggplot(res,aes(x=ID,y=enrichmentScore,fill=pvalue))+geom_col()+coord_flip()+
      scale_fill_gradientn(colors=c(color_low,color_mid,color_high)) + theme_classic() + theme(plot.margin = unit(c(2,2,2,2), "cm"))+
      theme(text=element_text(size=12,family="Arial"))
  } 
  
  if(type=="hyperG") {
    res <- res[res$enrich.ratio>absERcut&res$pvalue<pvalcut&res$p.adjust<fdrcut,]
    res <- res[order(res$enrich.ratio),]
    res$ID <- factor(res$Description,levels=res$Description,ordered = T)
    p <- ggplot(res,aes(x=ID,y=enrich.ratio,fill=pvalue))+geom_col()+coord_flip()+
      scale_fill_gradientn(colors=c(color_low,color_mid,color_high)) + theme_classic() + theme(plot.margin = unit(c(2,2,2,2), "cm"))+
      theme(text=element_text(size=12,family="Arial"))
  } 
  
  plotly_build(p)
}

dotplot <- function(res, absERcut=0.05,pvalcut=0.05,fdrcut=0.05,color_low="#51adcf",
                    color_mid="grey",color_high="#d3de32",type) {
  
  
  
  if(type=="gsea") {
    res <- res[res$pvalue<pvalcut&res$p.adjust<fdrcut,]
    res <- res[order(res$enrichmentScore),]
    res$ID <- factor(res$Description,levels=res$Description,ordered = T)
    p <- ggplot(res,aes(x=ID,y=enrichmentScore,color=pvalue))+geom_point()+coord_flip()+
      scale_color_gradientn(colors=c(color_low,color_mid,color_high)) + theme_classic() + theme(plot.margin = unit(c(2,2,2,2), "cm"))+
      theme(text=element_text(size=12,family="Arial"))
  } 
  
  if(type=="hyperG") {
    res <- res[res$enrich.ratio>absERcut&res$pvalue<pvalcut&res$p.adjust<fdrcut,]
    res <- res[order(res$enrich.ratio),]
    res$ID <- factor(res$Description,levels=res$Description,ordered = T)
    p <- ggplot(res,aes(x=ID,y=enrich.ratio,color=pvalue))+geom_point()+coord_flip()+
      scale_color_gradientn(colors=c(color_low,color_mid,color_high)) + theme_classic() + theme(plot.margin = unit(c(2,2,2,2), "cm"))+
      theme(text=element_text(size=12,family="Arial"))
  } 
  
  plotly_build(p)
}


network.plot <- function(res, absERcut=0.05,pvalcut=1,fdrcut=1,gene.list=NULL,color_low="#51adcf",
                    color_mid="grey",color_high="#d3de32",type) {

  if(type=="gsea") {
    geneList <- gene.list[,2]
    names(geneList) <- gene.list[,1]
    geneList <- sort(geneList,decreasing = T)
    y=res[res$pvalue < pvalcut,asis=T]
    y <- res[res$pvalue<pvalcut&res$p.adjust<fdrcut,asis=T]
    p <- cnetplot(y,colorEdge=T,foldChange = geneList)+scale_color_gradientn(colors=c(color_low,color_mid,color_high)) + theme_classic() + theme(plot.margin = unit(c(2,2,2,2), "cm"))+labs(x="",y="",color="fold change",size="counts")
  } 
  
  if(type=="hyperG") {
    res$ID <- res$Description;res$GeneRatio <- rep(1,nrow(res));res$BgRatio <- rep(1,nrow(res));
    res <- res[,c(9,1,10:11,2,7:8,6,5)]
    colnames(res) <- c("ID","Description","GeneRatio","BgRatio","pvalue","p.adjust", "qvalue", "geneID","Count")
    rownames(res) <- res$ID
    res$geneID <- as.character(res$geneID)
    geneSets <-lapply(res$geneID,function(x){strsplit(x,"/")[[1]]});names(geneSets) <- res$Description
    y <- new("enrichResult",
             readable = FALSE,
             result = res,
             pvalueCutoff = pvalcut,
             pAdjustMethod = "BH",
             qvalueCutoff = fdrcut,
             organism = "UNKNOWN",
             ontology = "UNKNOWN",
             gene = "UNKNOWN",
             keytype = "UNKNOWN",
             universe = "UNKNOWN",
             gene2Symbol = character(0),
             geneSets = geneSets)
    p <- cnetplot(y,colorEdge=T,foldChange = NULL)+scale_color_gradientn(colors=c(color_low,color_mid,color_high))  + theme_classic() + theme(plot.margin = unit(c(2,2,2,2), "cm"))+labs(x="",y="",color="fold change",size="counts")
  }
  
  #plotly_build(p)
  p
}

GSEA.plot <- function(gseaResult,geneSetID,title=" ",anno.cols="#d3de32",ifP=TRUE) {

    p <- gsea.xf(
      x=gseaResult,
      geneSetID,#富集的ID编号
      title=title,
      color = anno.cols,#GSEA线条颜色
      base_size = 11,
      rel_heights = c(1.5, 0.5, 1),
      pvalue_table = ifP, #是否添加 pvalue table
      ES_geom = "line" #running enrichment score用先还是用点ES_geom = "dot"
    )
    
    #p+theme_classic()
    p
}

decode.col.full <- function(){
  pals <- list(
    choices = list(
      Default = list("ggplot2" = hue_pal()(9)),
      Viridis = list(
        "viridis" = viridis_pal(option = "viridis")(10),
        "magma" = viridis_pal(option = "magma")(10),
        "inferno" = viridis_pal(option = "inferno")(10),
        "plasma" = viridis_pal(option = "plasma")(10),
        "cividis" = viridis_pal(option = "cividis")(10)
      ),
      Diverging = list(
        "BrBG" = brewer_pal(palette = "BrBG")(11), 
        "PiYG" = brewer_pal(palette = "PiYG")(11), 
        "PRGn" = brewer_pal(palette = "PRGn")(11), 
        "PuOr" = brewer_pal(palette = "PuOr")(11), 
        "RdBu" = brewer_pal(palette = "RdBu")(11), 
        "RdGy" = brewer_pal(palette = "RdGy")(11), 
        "RdYlBu" = brewer_pal(palette = "RdYlBu")(11), 
        "RdYlGn" = brewer_pal(palette = "RdYlGn")(11), 
        "Spectral" = brewer_pal(palette = "Spectral")(11)
      ), 
      Qualitative = list(
        "Accent" = brewer_pal(palette = "Accent")(8),
        "Dark2" = brewer_pal(palette = "Dark2")(8), 
        "Paired" = brewer_pal(palette = "Paired")(12), 
        "Pastel1" = brewer_pal(palette = "Pastel1")(9), 
        "Pastel2" = brewer_pal(palette = "Pastel2")(8), 
        "Set1" = brewer_pal(palette = "Set1")(8), 
        "Set2" = brewer_pal(palette = "Set2")(8), 
        "Set3" = brewer_pal(palette = "Set3")(12)
      ),
      Sequential = list(
        "Blues" = brewer_pal(palette = "Blues")(9),
        "BuGn" = brewer_pal(palette = "BuGn")(9),
        "BuPu" = brewer_pal(palette = "BuPu")(9), 
        "GnBu" = brewer_pal(palette = "GnBu")(9), 
        "Greens" = brewer_pal(palette = "Greens")(9), 
        "Greys" = brewer_pal(palette = "Greys")(9), 
        "Oranges" = brewer_pal(palette = "Oranges")(9), 
        "OrRd" = brewer_pal(palette = "OrRd")(9), 
        "PuBu" = brewer_pal(palette = "PuBu")(9), 
        "PuBuGn" = brewer_pal(palette = "PuBuGn")(9), 
        "PuRd" = brewer_pal(palette = "PuRd")(9), 
        "Purples" = brewer_pal(palette = "Purples")(9), 
        "RdPu" = brewer_pal(palette = "RdPu")(9), 
        "Reds" = brewer_pal(palette = "Reds")(9), 
        "YlGn" = brewer_pal(palette = "YlGn")(9), 
        "YlGnBu" = brewer_pal(palette = "YlGnBu")(9), 
        "YlOrBr" = brewer_pal(palette = "YlOrBr")(9), 
        "YlOrRd" = brewer_pal(palette = "YlOrRd")(9)
      )
    ), 
    textColor = c(
      rep(c("white", "black"), times = c(23, 18))
    )
  )
  return(pals)
}

decode.col <- function(){
  pals <- list(
    choices = list(
        "viridis" = viridis_pal(option = "viridis")(10),
        "magma" = viridis_pal(option = "magma")(10),
        "inferno" = viridis_pal(option = "inferno")(10),
        "plasma" = viridis_pal(option = "plasma")(10),
        "cividis" = viridis_pal(option = "cividis")(10),
        "BrBG" = brewer_pal(palette = "BrBG")(11), 
        "PiYG" = brewer_pal(palette = "PiYG")(11), 
        "PRGn" = brewer_pal(palette = "PRGn")(11), 
        "PuOr" = brewer_pal(palette = "PuOr")(11), 
        "RdBu" = brewer_pal(palette = "RdBu")(11), 
        "RdGy" = brewer_pal(palette = "RdGy")(11), 
        "RdYlBu" = brewer_pal(palette = "RdYlBu")(11), 
        "RdYlGn" = brewer_pal(palette = "RdYlGn")(11), 
        "Spectral" = brewer_pal(palette = "Spectral")(11),
        "Accent" = brewer_pal(palette = "Accent")(8),
        "Dark2" = brewer_pal(palette = "Dark2")(8), 
        "Paired" = brewer_pal(palette = "Paired")(12), 
        "Pastel1" = brewer_pal(palette = "Pastel1")(9), 
        "Pastel2" = brewer_pal(palette = "Pastel2")(8), 
        "Set1" = brewer_pal(palette = "Set1")(8), 
        "Set2" = brewer_pal(palette = "Set2")(8), 
        "Set3" = brewer_pal(palette = "Set3")(12),
        "Blues" = brewer_pal(palette = "Blues")(9),
        "BuGn" = brewer_pal(palette = "BuGn")(9),
        "BuPu" = brewer_pal(palette = "BuPu")(9), 
        "GnBu" = brewer_pal(palette = "GnBu")(9), 
        "Greens" = brewer_pal(palette = "Greens")(9), 
        "Greys" = brewer_pal(palette = "Greys")(9), 
        "Oranges" = brewer_pal(palette = "Oranges")(9), 
        "OrRd" = brewer_pal(palette = "OrRd")(9), 
        "PuBu" = brewer_pal(palette = "PuBu")(9), 
        "PuBuGn" = brewer_pal(palette = "PuBuGn")(9), 
        "PuRd" = brewer_pal(palette = "PuRd")(9), 
        "Purples" = brewer_pal(palette = "Purples")(9), 
        "RdPu" = brewer_pal(palette = "RdPu")(9), 
        "Reds" = brewer_pal(palette = "Reds")(9), 
        "YlGn" = brewer_pal(palette = "YlGn")(9), 
        "YlGnBu" = brewer_pal(palette = "YlGnBu")(9), 
        "YlOrBr" = brewer_pal(palette = "YlOrBr")(9), 
        "YlOrRd" = brewer_pal(palette = "YlOrRd")(9)
    ), 
    textColor = c(
      rep(c("white", "black"), times = c(23, 18))
    )
  )
  return(pals)
}


gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin=0
  df$ymax=0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

tableGrob2 <- function(d, p = NULL) {
  d <- d[order(rownames(d)),]
  tp <- tableGrob(d)
  if (is.null(p)) {
    return(tp)
  }
  pcol <- unique(ggplot_build(p)$data[[1]][["colour"]])
  j <- which(tp$layout$name == "rowhead-fg")
  
  for (i in seq_along(pcol)) {
    tp$grobs[j][[i+1]][["gp"]] = gpar(col = pcol[i])
  }
  return(tp)
}

gseaScores <- getFromNamespace("gseaScores", "DOSE")

gsea.xf <- function(x, geneSetID, title=" ", color, base_size = 11,
                    rel_heights=c(1.5, .5, 1), subplots = 1:3, pvalue_table = FALSE, ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_line(colour = "grey92"),
          panel.grid.minor = element_line(colour = "grey92"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description), size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description), size=1, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer +
    theme(legend.position = c(.8, .8), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          plot.margin=margin(t=.2, r = .2, b=0, l=.2, unit="cm"))
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          plot.margin = margin(t=-.1, b=0,unit="cm"),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (length(geneSetID) == 1) {

    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    
    col=c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[unique(inv)])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      alpha=.9,
      inherit.aes=FALSE)
  }
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0), color="grey")
  p.pos <- p.pos + ylab("Ranked list metric") +
    xlab("Rank in Ordered Dataset") +
    theme(plot.margin=margin(t = -.1, r = .2, b=.2, l=.2, unit="cm"))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    
    pd <- pd[,-1]
    pd <- round(pd, 4)
    
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  #plotlist[[1]]+plotlist[[2]]+plotlist[[3]]+ plot_layout(ncol = 1, heights = rel_heights)
  #return(plotlist)
  #plot_grid.xf(plotlist = plotlist, ncol=1, align="v", heights=rel_heights)
  plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
}

