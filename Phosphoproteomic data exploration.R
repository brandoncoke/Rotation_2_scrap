################################################################################
#All sites differentially affected by cyclin B depletion.
#Hégarat et al 10.15252/embj.2020104419
#Cyclin A triggers Mitosis either via the Greatwall kinase pathway or Cyclin B
#Table EV2
################################################################################
if (!requireNamespace("readxl", quietly = TRUE))
  install.packages("readxl")
library(readxl)
AllPhosposites<-read_xlsx("~\\Rotation 2\\Raw data\\Allphospho.xlsx",sheet=2)
Significantphospho<-read_xlsx("~\\Rotation 2\\Raw data\\Significantphospho.xlsx",sheet=2)
#Significantphospho= subset(Significantphospho,mean.fc > log10(2))
################################################################################
#GO analysis via gprofiler using gProfiler
################################################################################
if("gprofiler2" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("gprofiler2")}
library(gprofiler2)
gostres <- gost(query = unique(Significantphospho$Gene.names),
                organism = "hsapiens",sources="GO:BP",
                correction_method= "gSCS",significant=TRUE)
PhospoGO= gostres$result
publish_gosttable(gostres, highlight_terms = PhospoGO$term_id[1:15],
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", 
                                   "intersection_size"), filename = NULL)
################################################################################
#GO analysis via gprofiler using EnrichR. Package uses an API; queries 
################################################################################
if (!requireNamespace("enrichR", quietly = TRUE)){
  install.packages("enrichR")}
library(enrichR)
setEnrichrSite("Enrichr") # Human genes
dbs <- listEnrichrDbs()
GOterms <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "GO_Biological_Process_2018")
`%notin%` <- Negate(`%in%`)
genes= Significantphospho$Gene.names
ifelse(any(GOterms %notin% dbs$libraryName),
       print("INVALID LIBRARIES SELECTED"),
       enriched <- enrichr(genes, GOterms))
  
  
  
  
  
  
  
  
  
  
  
  
  