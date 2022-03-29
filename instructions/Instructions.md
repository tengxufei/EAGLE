---
output:
  html_document:
    theme: united
---

**The EAGLE shiny app allows users to analysis and visualize functional enrichment results, also provides gene annotation function.**

- **Upload** your data in the "Enrichment Analysis" tab and construct a functional enrichment analysis.
- **Visualize** the enrichment results in the "Enrichment Visualization" tab.
- **Annotate** your interested genes in the "Gene Annotation" tab.

<!-- 
-->
<a name="Instructions"></a> 

## Instructions

<!-- 
The app is hosted on the website:https://.shinyapps.io/EAGLE/ 
Code can be found on github: https://github.com/tengxufei/EAGLE
-->
Please post issues on [github](https://github.com/tengxufei/EAGLE/).

To run this app locally on your machine, download R or RStudio and run the following commands once to set up the environment:

```
install.packages(c("reshape2","ggplot2","ggthemes","gplots","dplyr","tidyr","DT","RColorBrewer","shinyBS","plotly","markdown","scales"))
```
You may now run the shiny app with just one command in R:
```
shiny::runGitHub("EAGLE", "tengxufei")
```

<a name="analysis"></a> 

## Enrichment Analysis

You may use this app to construct functional analysis by

1. Over-Representation Analysis
2. Gene Set Enrichment Analysis
3. Subpathway Analysis

<a name="input"></a>
### Input Data Format 

- For Over-Representation Analysis, users need to provide a list of their interested genes and a list of background genes, if the list of background genes is not provided, all annotated genes in selected functional datasets will be considered (We strongly suggest a upload list of background genes to avoid virtual low *p* value). Note that the gene lists must submitted in a gene symbol format (For human, HGNC symbols are suggested).

- For GSEA analysis, users need to provide a matrix contains their interested genes (first column) and their scores (universally the fold change of methylation or expression value in two or more conditions).

<a name="sets"></a>
### Functional data sets

#### Targets of m6A regulators
Potential targets of m6A regulators (writers, erasers, and readers) indicated by high-throughput sequencing such as CLIP-Seq, RIP-seq and ChIP-seq, etc.

#### Interaction targets of lncRNA
Potential targets of lncRNAs are
accessed from [RNAInter](http://www.rna-society.org/rnainter/).

#### Substrate of kinase
Substrate of kinase are retrieved from [PhosphoSitePlus](https://www.phosphosite.org/homeAction.action).

#### MSigDB (Molecular Signatures Database)
The description of MSigDB Collections is provided here: <http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp>



<a name="result"></a>
### Functional enrichment results

- For Over-Representation Analysis, the enrichment results contains 6 columns, which are "ID", "pvalue", "odds.ratio", "enrich.ratio", "counts", and "p.adjust".

- For GSEA analysis, the enrichment results contains 6 columns, which are "ID", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "rank", "leading_edge", and "core_enrichment".

<a name="visualization"></a>
## Visualization

### Bar Plot

This plot displays enriched functional items identified from hypergeometric test or GSEA analysis which match your cutoff.

### GSEA plot

This plot displays a standard GSEA plot for one or more enriched functional items identified from GSEA analysis which match your cutoff.

### Network Plot

This plot displays a item-gene network for enriched functional items identified from hypergeometric test or GSEA analysis which match your cutoff.



## Gene Annotation

### Based on functional gene sets

Users can upload a list of genes of interest (in m6A modification related research, usually modified genes or differentially modified genes) on the "Gene Annotation" page, and select "Functional Gene Sets" under the item "Annotation method". EAGLE will output an integration result including gene type, gene  description, possible upstream m6A regulators, KEGG pathway and MSigDB gene set to which the genes belong in tabular form.

### Based on literature mining

Users can select "Literature Mining" under the "Annotation method" item, and input the gene list of interest. EAGLE will give an integrated information related to users' input genes, including other genes, related chemicals,  diseases, and species. Those annotation data are retrieved from [PubTatorCentral](ftp://ftp.ncbi.nlm.nih.gov/pub/lu/PubTatorCentral/bioconcepts2pubtatorcentral.gz), only articles contains "DNA repair|m6A modification|Cell cycle|cell apoptosis|cell migration|cell proliferation|tumor development|cancer progression|cancer biomarker|tumorigenesis|P53|MYC" are included.  
