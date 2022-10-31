################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("BiocManager")}
library(Biobase)
if("dagLogo" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("dagLogo")}
library(dagLogo)
if("biomaRt" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("biomaRt")}
library(biomaRt)
if("UniProt.ws" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("UniProt.ws")}
library(UniProt.ws)
if("stringr" %in% rownames(installed.packages()) == FALSE){  
  install.packages("stringr")}
library(stringr)
################################################################################
#Data read and seperation into mitotic kinases
################################################################################
All_phosphosites= read.delim(
  "~/Rotation 2/Raw data/Kinase dataset from Phosphosite.tsv")
MitoticCDKsubstrates= subset(All_phosphosites,KINASE== "CDK1" |
                            KINASE== "CDK2")
Plk1substrates= All_phosphosites[All_phosphosites$KINASE== "PLK1",]
#Nek groupings based on 10.7554/eLife.44635.001
Group1NEKsubstrates= subset(All_phosphosites,KINASE== "NEK1" |
                              KINASE== "NEK3" |
                              KINASE== "NEK4")
#Group2NEKsubstrates= subset(All_phosphosites,KINASE== "NEK5" |
                              #KINASE== "NEK8")
Group3NEKsubstrates= subset(All_phosphosites,KINASE== "NEK2" |
                              KINASE== "NEK10")
Group4NEKsubstrates= subset(All_phosphosites,KINASE== "NEK1" |
                             KINASE== "NEK3" |
                              KINASE== "NEK4")
CK1substrates= All_phosphosites[grepl("CK1",All_phosphosites$KINASE),]
CK2substrates= All_phosphosites[grepl("CK2",All_phosphosites$KINASE),]
AKTsubstrates= All_phosphosites[grepl("Akt",All_phosphosites$KINASE),]
GSKsubstrates= All_phosphosites[grepl("GSK",All_phosphosites$KINASE),]
Aursubstrates= All_phosphosites[grepl("Aur",All_phosphosites$KINASE),]
rm(All_phosphosites)
################################################################################
#Creating extra column with AA sequence compatible for Daglogo
################################################################################
AA_for_Daglogo= function(string){
  first_halve= substr(string,1,8)
  first_halve= paste0(first_halve,"*")
  string= paste0(first_halve,substr(string,9,15))
  string= gsub('_','',string)
  string= toupper(string)
}
MitoticCDKsubstrates$Daglogo_input= lapply(MitoticCDKsubstrates$SITE_...7_AA, AA_for_Daglogo)
Plk1substrates$Daglogo_input= lapply(Plk1substrates$SITE_...7_AA, AA_for_Daglogo)
Group1NEKsubstrates$Daglogo_input= lapply(Group1NEKsubstrates$SITE_...7_AA, AA_for_Daglogo)
Group3NEKsubstrates$Daglogo_input= lapply(Group3NEKsubstrates$SITE_...7_AA, AA_for_Daglogo)
Group4NEKsubstrates$Daglogo_input= lapply(Group4NEKsubstrates$SITE_...7_AA, AA_for_Daglogo)
CK1substrates$Daglogo_input= lapply(CK1substrates$SITE_...7_AA, AA_for_Daglogo)
CK2substrates$Daglogo_input= lapply(CK2substrates$SITE_...7_AA, AA_for_Daglogo)
AKTsubstrates$Daglogo_input= lapply(AKTsubstrates$SITE_...7_AA, AA_for_Daglogo)
GSKsubstrates$Daglogo_input= lapply(GSKsubstrates$SITE_...7_AA, AA_for_Daglogo)
Aursubstrates$Daglogo_input= lapply(Aursubstrates$SITE_...7_AA, AA_for_Daglogo)
################################################################################
#Getting z-scores for amino acids adjacent to sites
################################################################################
mart <- useMart("ensembl")
humanmart <-
  useDataset(mart = mart, dataset = "hsapiens_gene_ensembl")
Getphospeps=function(ID,anchorPos,range=7){
  output <- fetchSequence(IDs = IDs,
                       anchorAA = "*",
                       anchorPos = anchorPos,
                       mart = humanmart,
                       type="hgnc_symbol",
                       upstreamOffset = range,
                       downstreamOffset = range)
  
  return(output)
}
anchorPos=as.character(MitoticCDKsubstrates$Daglogo_input)
IDs=as.character(MitoticCDKsubstrates$SUBSTRATE)
seqs= Getphospeps(IDs,
                  anchorPos,
                  range=7)
proteome <- dagLogo::prepareProteome("UniProt", species = "Homo sapiens")
bg_ztest <- buildBackgroundModel(seqs, background = "wholeProteome", 
                                 proteome = proteome, testType = "ztest")
t0 <- testDAU(seqs, dagBackground = bg_ztest)
Substrate=c(MitoticCDKsubstrates$SUBSTRATE,CK1substrates$SUBSTRATE,
            CK2substrates$SUBSTRATE,Group1NEKsubstrates$SUBSTRATE,
            GSKsubstrates$SUBSTRATE,AKTsubstrates$SUBSTRATE,
            Plk1substrates$SUBSTRATE)
Kinase=c(MitoticCDKsubstrates$KINASE,
         CK1substrates$KINASE,
         CK2substrates$KINASE,Group1NEKsubstrates$KINASE,
         GSKsubstrates$KINASE,
         AKTsubstrates$KINASE,
         Plk1substrates$KINASE)
NNdataframe=data.frame(AA1_zscore="X",
                       AA2_zscore="X",
                       AA3_zscore="X",
                       AA4_zscore="X",
                       AA5_zscore="X",
                       AA6_zscore="X",
                       AA7_zscore="X",
                       AA8_zscore="X",
                       AA9_zscore="X",
                       AA10_zscore="X",
                       AA11_zscore="X",
                       AA12_zscore="X",
                       AA13_zscore="X",
                       AA14_zscore="X",
                       Substrate,
                       Kinase)
sequences=c(as.character(MitoticCDKsubstrates$Daglogo_input),
            as.character(CK1substrates$Daglogo_input),
            as.character(CK2substrates$Daglogo_input),
            as.character(Group1NEKsubstrates$Daglogo_input),
            as.character(GSKsubstrates$Daglogo_input),
            as.character(AKTsubstrates$Daglogo_input),
            as.character(Plk1substrates$Daglogo_input))
sequences= gsub('[[:punct:] ]+', 'X', sequences)

c=1
cycles=length(sequences) + 1
while(c<cycles){
one_up_from_target=str_locate(sequences[c],"X")
one_up_from_target=one_up_from_target[1]
NNdataframe$AA7_zscore[c]=substr(sequences[c],one_up_from_target-1,one_up_from_target-1)
NNdataframe$AA6_zscore[c]=substr(sequences[c],one_up_from_target-2,one_up_from_target-2)
NNdataframe$AA5_zscore[c]=substr(sequences[c],one_up_from_target-3,one_up_from_target-3)
NNdataframe$AA4_zscore[c]=substr(sequences[c],one_up_from_target-4,one_up_from_target-4)
NNdataframe$AA3_zscore[c]=substr(sequences[c],one_up_from_target-5,one_up_from_target-5)
NNdataframe$AA2_zscore[c]=substr(sequences[c],one_up_from_target-6,one_up_from_target-6)
NNdataframe$AA1_zscore[c]=substr(sequences[c],one_up_from_target-7,one_up_from_target-7)
    
NNdataframe$AA8_zscore[c]=substr(sequences[c],one_up_from_target+1,one_up_from_target+1)
NNdataframe$AA9_zscore[c]=substr(sequences[c],one_up_from_target+2,one_up_from_target+2)
NNdataframe$AA10_zscore[c]=substr(sequences[c],one_up_from_target+3,one_up_from_target+3)
NNdataframe$AA11_zscore[c]=substr(sequences[c],one_up_from_target+4,one_up_from_target+4)
NNdataframe$AA12_zscore[c]=substr(sequences[c],one_up_from_target+5,one_up_from_target+5)
NNdataframe$AA13_zscore[c]=substr(sequences[c],one_up_from_target+6,one_up_from_target+6)
NNdataframe$AA14_zscore[c]=substr(sequences[c],one_up_from_target+7,one_up_from_target+7)
NNdataframe
c=c+1}
zscore_matrix= t0@statistics
zscore_matrix= rbind(zscore_matrix,X=0)
letters=rownames(zscore_matrix)
NNdataframe$AA1_zscore[NNdataframe$AA1_zscore == ""] = "X"
NNdataframe$AA2_zscore[NNdataframe$AA2_zscore == ""] = "X"
NNdataframe$AA3_zscore[NNdataframe$AA3_zscore == ""] = "X"
NNdataframe$AA4_zscore[NNdataframe$AA4_zscore == ""] = "X"
NNdataframe$AA5_zscore[NNdataframe$AA5_zscore == ""] = "X"
NNdataframe$AA6_zscore[NNdataframe$AA6_zscore == ""] = "X"
NNdataframe$AA7_zscore[NNdataframe$AA7_zscore == ""] = "X"
NNdataframe$AA8_zscore[NNdataframe$AA8_zscore == ""] = "X"
NNdataframe$AA9_zscore[NNdataframe$AA9_zscore == ""] = "X"
NNdataframe$AA10_zscore[NNdataframe$AA10_zscore == ""] = "X"
NNdataframe$AA11_zscore[NNdataframe$AA11_zscore == ""] = "X"
NNdataframe$AA12_zscore[NNdataframe$AA12_zscore == ""] = "X"
NNdataframe$AA13_zscore[NNdataframe$AA13_zscore == ""] = "X"
NNdataframe$AA14_zscore[NNdataframe$AA14_zscore == ""] = "X"

c=1
dict= t(zscore_matrix[,8])
while(c<cycles){
  index=which(NNdataframe[c,7] == letters)
  NNdataframe[c,7]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,9])
while(c<cycles){
  index=which(NNdataframe[c,8] == letters)
  NNdataframe[c,8]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,10])
while(c<cycles){
  index=which(NNdataframe[c,9] == letters)
  NNdataframe[c,9]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,11])
while(c<cycles){
  index=which(NNdataframe[c,10] == letters)
  NNdataframe[c,10]= dict[index]
  c=c+1
}



c=1
dict= t(zscore_matrix[,12])
while(c<cycles){
  index=which(NNdataframe[c,11] == letters)
  NNdataframe[c,11]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,13])
while(c<cycles){
  index=which(NNdataframe[c,12] == letters)
  NNdataframe[c,12]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,14])
while(c<cycles){
  index=which(NNdataframe[c,13] == letters)
  NNdataframe[c,13]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,15])
while(c<cycles){
  index=which(NNdataframe[c,14] == letters)
  NNdataframe[c,14]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,7])
while(c<cycles){
  index=which(NNdataframe[c,6] == letters)
  NNdataframe[c,6]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,6])
while(c<cycles){
  index=which(NNdataframe[c,5] == letters)
  NNdataframe[c,5]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,5])
while(c<cycles){
  index=which(NNdataframe[c,4] == letters)
  NNdataframe[c,4]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,4])
while(c<cycles){
  index=which(NNdataframe[c,3] == letters)
  NNdataframe[c,3]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,3])
while(c<cycles){
  index=which(NNdataframe[c,2] == letters)
  NNdataframe[c,2]= dict[index]
  c=c+1
}

c=1
dict= t(zscore_matrix[,2])
while(c<cycles){
  index=which(NNdataframe[c,1] == letters)
  NNdataframe[c,1]= dict[index]
  c=c+1
}

write.csv(NNdataframe,"K:\\CDK1 z-score matrix.csv")
