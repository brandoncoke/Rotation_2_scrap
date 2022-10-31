################################################################################
#Package install
################################################################################
if("VennDiagram" %in% rownames(installed.packages()) == FALSE){
  install.packages("VennDiagram")}
library(VennDiagram)
if("stringr" %in% rownames(installed.packages()) == FALSE){
  install.packages("stringr")}
library(stringr)
if("readxl" %in% rownames(installed.packages()) == FALSE){
  install.packages("readxl")}
library(readxl)
################################################################################
#Load Sequences of interest.RDS...
################################################################################
Known_phosphosites <- read.delim("~/Rotation 2/Raw data/Kinase dataset from Phosphosite.tsv")
Phospho_sequences= read_xlsx("~\\Rotation 2\\Raw data\\Phosphosites significantly changed by cyclin B depeletion.xlsx",sheet=2)
Known_phosphosites= subset(Known_phosphosites,SUBSTRATE %in% Phospho_sequences$Gene.names)
Known_phosphosites$SUB_MOD_RSD= gsub("[^0-9.-]", "", Known_phosphosites$SUB_MOD_RSD)
Phospho_sequences$Already_known= F
Phospho_sequences$Known_kinase= NA
c=1
cycles=nrow(Known_phosphosites) + 1
while(c<cycles){
  target= Known_phosphosites$SUBSTRATE[c]
  site= Known_phosphosites$SUB_MOD_RSD[c]
  tempframe= subset(Phospho_sequences,Gene.names == target)
  if(any(grepl(site,tempframe$Position.in.peptide))){
    site_to_be_changed= tempframe$phosphos.id[grepl(site,tempframe$Position.in.peptide)]
    site_to_be_changed=site_to_be_changed[1]
    Phospho_sequences$Already_known[
      Phospho_sequences$phosphos.id == site_to_be_changed] = T
    Phospho_sequences$Known_kinase[
      Phospho_sequences$phosphos.id == site_to_be_changed] = Known_phosphosites$KINASE[c]
  }
  
  c=c+1
}