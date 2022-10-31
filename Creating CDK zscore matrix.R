################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("Biobase")}
library(BiocManager)
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
#Read of raw phosphosite data
################################################################################
All_phosphosites= read.delim(
  "~/Rotation 2/Raw data/Kinase dataset from Phosphosite.tsv")
################################################################################
#Creating basic zscore matrix
################################################################################
CDK1only= All_phosphosites[All_phosphosites$KINASE == "CDK1",]
CDKonly= All_phosphosites[grepl("CDK",All_phosphosites$KINASE),]
AA_matrix=data.frame(AA1="X",
                       AA2="X",
                       AA3="X",
                       AA4="X",
                       AA5="X",
                       AA6="X",
                       AA7="X",
                       AA8="X",
                       AA9="X",
                       AA10="X",
                       AA11="X",
                       AA12="X",
                       AA13="X",
                       AA14="X",
                       Substrate=CDK1only$SUBSTRATE,
                       Kinase= CDK1only$KINASE,
                       Organism=CDK1only$KIN_ORGANISM)
################################################################################
#Marking +1 of phosphosequences with a X using AA_for_matrix function
#then inserting phosphosequences into zscore matrix
################################################################################
AA_for_matrix= function(string){
  first_halve= substr(string,1,8)
  first_halve= paste0(first_halve,"X")
  string= paste0(first_halve,substr(string,9,15))
  string= gsub('_','',string)
  string= toupper(string)
  string
}

sequences= AA_for_matrix(CDK1only$SITE_...7_AA)
for(i in 1:length(sequences)){
  one_up_from_target=str_locate(sequences[i],"X")
  one_up_from_target=one_up_from_target[1]
  AA_matrix$AA7[i]=substr(sequences[i],one_up_from_target-1,one_up_from_target-1)
  AA_matrix$AA6[i]=substr(sequences[i],one_up_from_target-2,one_up_from_target-2)
  AA_matrix$AA5[i]=substr(sequences[i],one_up_from_target-3,one_up_from_target-3)
  AA_matrix$AA4[i]=substr(sequences[i],one_up_from_target-4,one_up_from_target-4)
  AA_matrix$AA3[i]=substr(sequences[i],one_up_from_target-5,one_up_from_target-5)
  AA_matrix$AA2[i]=substr(sequences[i],one_up_from_target-6,one_up_from_target-6)
  AA_matrix$AA1[i]=substr(sequences[i],one_up_from_target-7,one_up_from_target-7)
  
  AA_matrix$AA8[i]=substr(sequences[i],one_up_from_target+1,one_up_from_target+1)
  AA_matrix$AA9[i]=substr(sequences[i],one_up_from_target+2,one_up_from_target+2)
  AA_matrix$AA10[i]=substr(sequences[i],one_up_from_target+3,one_up_from_target+3)
  AA_matrix$AA11[i]=substr(sequences[i],one_up_from_target+4,one_up_from_target+4)
  AA_matrix$AA12[i]=substr(sequences[i],one_up_from_target+5,one_up_from_target+5)
  AA_matrix$AA13[i]=substr(sequences[i],one_up_from_target+6,one_up_from_target+6)
  AA_matrix$AA14[i]=substr(sequences[i],one_up_from_target+7,one_up_from_target+7)}
################################################################################
#Ensuring empty positions contain X
################################################################################
AA_matrix$AA1[AA_matrix$AA1 == ""] = "X"
AA_matrix$AA2[AA_matrix$AA2 == ""] = "X"
AA_matrix$AA3[AA_matrix$AA3 == ""] = "X"
AA_matrix$AA4[AA_matrix$AA4 == ""] = "X"
AA_matrix$AA5[AA_matrix$AA5 == ""] = "X"
AA_matrix$AA6[AA_matrix$AA6 == ""] = "X"
AA_matrix$AA7[AA_matrix$AA7 == ""] = "X"
AA_matrix$AA8[AA_matrix$AA8 == ""] = "X"
AA_matrix$AA9[AA_matrix$AA9 == ""] = "X"
AA_matrix$AA10[AA_matrix$AA10 == ""] = "X"
AA_matrix$AA11[AA_matrix$AA11 == ""] = "X"
AA_matrix$AA12[AA_matrix$AA12 == ""] = "X"
AA_matrix$AA13[AA_matrix$AA13 == ""] = "X"
AA_matrix$AA14[AA_matrix$AA14 == ""] = "X"
################################################################################
#Removing any rows with invalid residues in 7 position i.e. phoshporylated
#residue not T, S or Y.
################################################################################
valid_phosphorylated_residues= c("T","S","Y")
AA_matrix= AA_matrix[
  AA_matrix$AA7 %in% valid_phosphorylated_residues,]

rm(valid_phosphorylated_residues)
################################################################################
#Calculating z-scores for the residues at a particular position
#defining zscore calculation first
################################################################################
zscore= function(x,a_vector){ #Calculates standard zscore
  mean=mean(a_vector)
  sd= sd(a_vector)
  (x-mean)/sd
}
`%notin%` <- Negate(`%in%`)
AAs=c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S",
             "T","V","W","X","Y")
obtain_zscores= function(AA_letters){
  df= as.data.frame(table(AA_letters)) #converts table to dataframe
  AA_vector= c(t(df[2])) #Getting the frequency of the AA in same order specificed by AAs
  df$zscore= zscore(AA_vector,AA_vector) #Obtains the zscore of the amino acids 
  if(length(AA_letters[AAs %notin% df$AA_letters]) != 0){ #Redundant check performed previously when debugging
    fill_data_frame=data.frame(AA_letters= AAs[AAs %notin% df$AA_letters], Freq= NA,
                               zscore=NA)
    df= rbind(df,fill_data_frame)}
  df[,-2]
}

zscore_matrix= data.frame(AA_letters= AAs)
################################################################################
#Obtaining the zscores for each position
################################################################################
AA1_zscores= obtain_zscores(AA_matrix$AA1)
AA2_zscores= obtain_zscores(AA_matrix$AA2)
AA3_zscores= obtain_zscores(AA_matrix$AA3)
AA4_zscores= obtain_zscores(AA_matrix$AA4)
AA5_zscores= obtain_zscores(AA_matrix$AA5)
AA6_zscores= obtain_zscores(AA_matrix$AA6)
AA7_zscores= obtain_zscores(AA_matrix$AA7)
AA8_zscores= obtain_zscores(AA_matrix$AA8)
AA9_zscores= obtain_zscores(AA_matrix$AA9)
AA10_zscores= obtain_zscores(AA_matrix$AA10)
AA11_zscores= obtain_zscores(AA_matrix$AA11)
AA12_zscores= obtain_zscores(AA_matrix$AA12)
AA13_zscores= obtain_zscores(AA_matrix$AA13)
AA14_zscores= obtain_zscores(AA_matrix$AA14)
################################################################################
#Merging the zscores to zscore matrix
################################################################################
zscore_matrix= merge(zscore_matrix,AA1_zscores)
colnames(zscore_matrix)[2]= "AA1_zscore"
zscore_matrix= merge(zscore_matrix,AA2_zscores)
colnames(zscore_matrix)[3]= "AA2_zscore"
zscore_matrix= merge(zscore_matrix,AA3_zscores)
colnames(zscore_matrix)[4]= "AA3_zscore"
zscore_matrix= merge(zscore_matrix,AA4_zscores)
colnames(zscore_matrix)[5]= "AA4_zscore"
zscore_matrix= merge(zscore_matrix,AA5_zscores)
colnames(zscore_matrix)[6]= "AA5_zscore"
zscore_matrix= merge(zscore_matrix,AA6_zscores)
colnames(zscore_matrix)[7]= "AA6_zscore"
zscore_matrix= merge(zscore_matrix,AA7_zscores)
colnames(zscore_matrix)[8]= "AA7_zscore"
zscore_matrix= merge(zscore_matrix,AA8_zscores)
colnames(zscore_matrix)[9]= "AA8_zscore"
zscore_matrix= merge(zscore_matrix,AA9_zscores)
colnames(zscore_matrix)[10]= "AA9_zscore"
zscore_matrix= merge(zscore_matrix,AA10_zscores)
colnames(zscore_matrix)[11]= "AA10_zscore"
zscore_matrix= merge(zscore_matrix,AA11_zscores)
colnames(zscore_matrix)[12]= "AA11_zscore"
zscore_matrix= merge(zscore_matrix,AA12_zscores)
colnames(zscore_matrix)[13]= "AA12_zscore"
zscore_matrix= merge(zscore_matrix,AA13_zscores)
colnames(zscore_matrix)[14]= "AA13_zscore"
zscore_matrix= merge(zscore_matrix,AA14_zscores)
colnames(zscore_matrix)[15]= "AA14_zscore"
################################################################################
#Exporting
################################################################################
write.csv(zscore_matrix,"~\\CDK1 matrix.csv")
