################################################################################
#Package install
################################################################################
if("readxl" %in% rownames(installed.packages()) == FALSE){
  install.packages("readxl")}
library(readxl)
################################################################################
#loading and removing excess phosphoproteomic data i.e. proteins not of interest.
################################################################################
POIs= c("FAM207A",
  "FAM208B",
  "HN1L",
  "KIAA1143",
  "KIAA1671",
  "PDXDC1",
  "PRRC2C",
  "RBMS2") #POIs proteins of interest
sequences_and_genes <- read_xlsx(
  "~\\Rotation 2\\Raw data\\All of the phosphosites\\All phosphorylation sites affected by cyclin B depeletion.xlsx",sheet=2)
Other_POIs <- read_xlsx(
  "~\\Rotation 2\\Raw data\\Proteins of interest\\Raw POIs.xlsx",sheet=3)
POIs= subset(sequences_and_genes,Gene.names %in% POIs)
write.csv(POIs,"~\\POIs.csv")
Other_POIs=subset(sequences_and_genes,Gene.names %in% Other_POIs$`Gene name`)
write.csv(Other_POIs,"~\\Other_POIs.csv")
