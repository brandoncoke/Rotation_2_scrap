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
if("stringr" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("stringr")}
library(stringr)
################################################################################
#Loading raw data of POIs and other proteins affected by cyclin B depeletion
################################################################################
POIs <- read.csv("~/Rotation 2/Raw data/Proteins of interest fastas and alignments/POIs.csv")
Other_POIs <- read.csv("~/Rotation 2/Raw data/Proteins of interest fastas and alignments/Other_POIs.csv")
POIs$Gene.names[POIs$Gene.names == "HN1L"] = "JPT2"
POIs$focused= T
Other_POIs$focused= F
POIs= rbind(POIs,Other_POIs)
rm(Other_POIs)
POIs$Gene.names[POIs$Gene.names == "HN1L"] = "JPT2" #As protein as another name
################################################################################
#Getting full protein sequences usng biomaRt
################################################################################
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
IDs= gsub(";.*","",POIs$Gene.names)
IDs = as.character(unique(IDs))
Sequences = getSequence(id=as.character(IDs),
                        type="hgnc_symbol",
                        seqType="peptide", 
                        mart=mart)
rm(IDs,mart)
Sequences$peptide=gsub("X",NA,Sequences$peptide)
Sequences$peptide=gsub("Sequence unavailable",NA,Sequences$peptide)
Sequences$Peptide_length=nchar(Sequences$peptide)
Sequences=na.omit(Sequences)
#library(dplyr)
#All_peptide_sequences %>% group_by(Peptide_length) %>% top_n(1, Peptide_length)
c=2
genes=unique(Sequences$hgnc_symbol)
tempframe= subset(Sequences,hgnc_symbol ==genes[1])
Cleaned_sequences=tempframe[tempframe$Peptide_length == max(tempframe$Peptide_length),]
cycles=length(genes) +1
while(c<cycles){
  tempframe= subset(Sequences,hgnc_symbol ==genes[c])
  tempframe=tempframe[tempframe$Peptide_length == max(tempframe$Peptide_length),]
  tempframe=tempframe[1,]
  Cleaned_sequences=rbind(Cleaned_sequences,tempframe)
  tempframe= NA
  c=c+1
}
rm(tempframe,genes,Sequences)
################################################################################
#Cleaning phosphopeptide sequences for analysis i.e. removing probability brackets
################################################################################
Booleans_of_POIs= POIs[,30:50]
Fold_changes= POIs$mean.fc
p_value= POIs$p.val
#POIs$entrezgene_id= as.factor(All_peptide_POIs$hgnc_symbol)
POIs$raw_phospho= POIs$Phospho..STY..Probabilities
POIs$Phospho..STY..Probabilities= 
  gsub('[[:digit:]]+', '', POIs$Phospho..STY..Probabilities)
POIs$Phospho..STY..Probabilities= 
  gsub('[[:punct:] ]+', '', POIs$Phospho..STY..Probabilities)
colnames(POIs)[6]= "Phosphorylation_site"
POIs$Protein.names= gsub(";.*","",POIs$Gene.names)
POIs= data.frame(Gene_name=POIs$Gene.names,
                              Phospo_POIs=POIs$Phosphorylation_site,
                 Raw_phosphoso=POIs$raw_phospho,
                              Surrounding_AA=NA)
POIs= cbind(POIs,Fold_changes,p_value,Booleans_of_POIs)
rm(Booleans_of_POIs)
POIs$Potential_Polobox=grepl(
  "STS",POIs$Phospo_POIs) #Identifying potential polobox binding motifs
`%notin%` <- Negate(`%in%`)
POIsnosequence=subset(POIs,Gene_name %notin% Cleaned_sequences$hgnc_symbol)
POIs=subset(POIs,Gene_name %in% Cleaned_sequences$hgnc_symbol)
################################################################################
#Adding full protein sequences to POIs dataframe
################################################################################
c=1
cycles=nrow(POIs) +1
while(c<cycles){
  id= POIs$Gene_name[c]
  POIs$Full_sequence[c]=
    Cleaned_sequences$peptide[Cleaned_sequences$hgnc_symbol == id]
  c=c+1
}
################################################################################
#Obtaining residues within 35 either side of the phosphopeptide.
################################################################################
c=1
while(c<cycles){ 
  full_sequence=POIs$Full_sequence[c]
  Phosphosequence=POIs$Phospo_POIs[c]
  residue_location= stringr::str_locate(full_sequence,Phosphosequence)
  residue_location=residue_location[1]
  new_sequence= substr(full_sequence,residue_location-35,residue_location+35)
  POIs$Surrounding_AA[c]=new_sequence
  c=c+1}
POIsnosequence$Full_sequence= NA
POIs= rbind(POIs,POIsnosequence)
rm(c,cycles,full_sequence,id,new_sequence,residue_location,
   Cleaned_sequences,POIsnosequence,Phosphosequence)

################################################################################
#Identfying cyclin binding motifs in the full protein sequences...
################################################################################
#Change to Surrounding_AA to assess if motif nearby
#Bascically uses grep to see if the sequences have a particular motif
################################################################################
#Cy motif
################################################################################
POIs$Contains_Cy_CDK_motif=FALSE
POIs$Contains_Cy_CDK_motif[grepl("R.L",POIs$Full_sequence)]= #As | aint working in regular expression
  TRUE
################################################################################
#G2 specific cyclin binding motif. Not valid only a yeast cyclin binding motif. 
################################################################################
POIs$Contains_PXXPXF_motif=FALSE
POIs$Contains_PXXPXF_motif[grepl("P..P.F",POIs$Full_sequence)]=
  TRUE
################################################################################
#M phase specific cyclin binding motif 10.1016/j.molcel.2019.04.026
################################################################################
POIs$Contains_LXF_motif=FALSE
POIs$Contains_LXF_motif[grepl("L.F",POIs$Full_sequence)]= 
  TRUE
################################################################################
#CKS binding motif
################################################################################
POIs$Contains_CKS_binding_consensus=F
POIs$Contains_CKS_binding_consensus[grepl("F.TP",POIs$Full_sequence)]= 
  TRUE
POIs$Contains_CKS_binding_consensus[grepl("I.TP",POIs$Full_sequence)]=
  TRUE
POIs$Contains_CKS_binding_consensus[grepl("L.TP",POIs$Full_sequence)]= 
  TRUE
POIs$Contains_CKS_binding_consensus[grepl("P.TP",POIs$Full_sequence)]=
  TRUE
POIs$Contains_CKS_binding_consensus[grepl("V.TP",POIs$Full_sequence)]= 
  TRUE
POIs$Contains_CKS_binding_consensus[grepl("W.TP",POIs$Full_sequence)]=
  TRUE
################################################################################
#Last of cyclin binding motifs booleans 10.15252/embj.2020105839 and 10.1016/j.celrep.2020.107757
################################################################################
POIs$Contains_NLXXL= F
POIs$Contains_NLXXL[grepl("NL..L",POIs$Full_sequence)]=
  TRUE
POIs$Contains_S_phase= F
POIs$Contains_S_phase[grepl("(R|K).L(F|L)",POIs$Full_sequence)]=
  TRUE
POIs$Contains_Common= F
POIs$Contains_Common[grepl("(R|K).L.(F|L)",POIs$Full_sequence)]=
  TRUE
Focused_POIs= subset(POIs,focused == T)
################################################################################
#Export
################################################################################
write.csv(Focused_POIs,"~\\Processed POIs sequences + attributes.csv")
