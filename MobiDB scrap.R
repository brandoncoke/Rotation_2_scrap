################################################################################
#Package install 
################################################################################
if("jsonlite" %in% rownames(installed.packages()) == FALSE){
  install.packages("jsonlite")}
library(jsonlite)
if("httr" %in% rownames(installed.packages()) == FALSE){
  install.packages("httr")}
library(httr)

accessioncodes<- "acc=P10912,P04050"
endpoint <- "https://mobidb.bio.unipd.it/api/download?"
endpoint2 <- "&projection=acc,name,organism,ncbi_taxon_id,proteome,sequence,prediction-disorder-mobidb_lite.regions&format=json"
httr::GET(paste0(endpoint,accessioncodes,endpoint2),
          httr::timeout(10)) -> response
if (response$status_code >= 400){
  err_msg = httr::http_status(response)
  stop(err_msg); 
  paste("Error code: ", "response$status_code")}
content <- httr::content(response,as="text",encoding="UTF-8")
#content<- jsonlite::validate(content)

data = jsonlite::fromJSON(content)