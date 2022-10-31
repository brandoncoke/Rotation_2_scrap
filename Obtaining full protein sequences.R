################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("BiocManager")}
library(Biobase)
if("readxl" %in% rownames(installed.packages()) == FALSE){
  install.packages("readxl")}
library(readxl)
if("biomaRt" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("biomaRt")}
library(biomaRt)
################################################################################
#Package install
################################################################################
sequences_and_genes <- read_xlsx("~\\Rotation 2\\Raw data\\Phosphosites significantly changed by cyclin B depeletion.xlsx",sheet=2)
#ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
IDs= gsub(";.*","",sequences_and_genes$Gene.names)
IDs = as.character(unique(IDs))
Sequences = getSequence(id=as.character(IDs),
                        type="hgnc_symbol",
                        seqType="peptide", 
                        mart=mart)
rm(ensembl,sequences_and_genes,Sequences_temp,IDs,mart)

save.image("~/Sequences of interest.RData") #Saving progress so far into working directory.
write.csv(Sequences,"K:\\All sequences.csv")

