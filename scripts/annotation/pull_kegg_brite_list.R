library(KEGGREST)
library(jsonlite)
library(tidyverse)

retrieveBRITE <- function(id,organismCode=""){
  query <- KEGGREST::keggGet(paste0("br:",id,"/json")) %>% fromJSON
  return(flattenJSON(query,organismCode))
}

flattenJSON <- function(json,organismCode="") {
  df <- convertChildren(json,parent="",organismCode) %>% bind_rows
  max_parents <- max(lengths(str_split(df$parent,pattern="_")))
  return(df%>%separate(parent,into=paste0("parent",seq(0,max_parents-1)),fill="right",
                       sep="_",remove=FALSE) %>%
           mutate(BRITE=parent1) %>%
           select(-parent0,-parent1))
}

convertChildren <- function(input,parent,organismCode,genesOut=list()) {
  if(is.null(input)) {return(genesOut)}
  t <- map_if(input,is.data.frame,list) %>% as_tibble
  if ("children" %in% colnames(t)) {
    for (i in 1:length(t$name)) {
      levelName <- paste0(parent,"_",t$name[[i]])
      print(levelName)
      genesOut <- convertChildren(t$children[[i]],levelName,organismCode,genesOut)
    }
  }
  else {
    for (i in 1:length(t$name)) {
      n <- str_split_1(t$name[i],pattern="\t")
      ids <- str_split_1(str_split_1(n[1],pattern=";")[1],pattern=" ")
      keggName <- paste0(organismCode,":",ids[1])
      geneName <- ids[2]
      ncbiName <- str_split_1(KEGGREST::keggConv("ncbi-proteinid",kegg_name),":")[2]
      genesOut <- append(genesOut,list(tibble("parent"=parent,
                                              "kegg_name"=keggName,
                                              "gene_name"=geneName,
                                              "ncbi_name"=ncbiName)))
    }
  }
  return(genesOut)
}
