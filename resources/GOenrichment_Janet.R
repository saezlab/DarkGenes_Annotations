library(Biobase)
# source("https://bioconductor.org/biocLite.R")
# biocLite("GOstats")
library(GOstats)
library(org.Hs.eg.db) ## annotation Package for Homo Sapiens 
library(annotate)
library(AnnotationDbi)
library(BiocGenerics)

haveGO = function(entrezIds) {
  allGO = org.Hs.egGO
  all_genes = mappedkeys(allGO)
  have = sapply(entrezIds, function(x){
    if (length(x) == 1 && is.na(x)){
      FALSE
    }
    else if (x%in% all_genes){
      TRUE
    }
    else{
      FALSE
    }
  }
  )
  filteredIds <- entrezIds[have]
}

filterbyGO = function(entrezIds, ontologyType=c("BP", "MF", "CC")) {
  print("filter by GO")
  chip = "org.Hs.eg.db"
  ontologyType=match.arg(ontologyType)
  
  haveGo <- sapply(mget(entrezIds, getAnnMap(map="GO", chip=chip, type="db")),
                   function(x) {
                     if (length(x) == 1 && is.na(x)){
                       FALSE
                     }
                     else {
                       onts <- subListExtract(x, "Ontology",
                                              simplify=TRUE)
                       ontologyType %in% onts
                     }
                     
                   })
  
  filteredIds <- names(haveGo)[haveGo]
  return(filteredIds)
}

perform_GO_Enrichment =  function(entrezIds, ontologyType=c("BP", "MF", "CC")) {
  
  # filter for GO branch
  genesBranch = filterbyGO(entrezIds, ontologyType )
  print("genes filterbyGO branch has length of:")
  length(genesBranch)
  
  if (ontologyType == "BP"){
    universe = universeBP
  }
  if (ontologyType == "MF"){
    universe = universeMF
  }
  
  # GO biological process
  params = new("GOHyperGParams", 
               geneIds=genesBranch, 
               universeGeneIds=universe,
               annotation="org.Hs.eg.db", 
               ontology=ontologyType , 
               pvalueCutoff=0.05, 
               conditional=TRUE, 
               testDirection="over")
  
  hgOver=hyperGTest(params) ## hyperGTest from GOstats
  sumTest <-summary(hgOver)
  sumTest$padj <- p.adjust(sumTest$Pvalue)
  return(sumTest)
}

