################################################################################
#Package install
################################################################################
if("stringr" %in% rownames(installed.packages()) == FALSE){  
  install.packages("stringr")}
library(stringr)
################################################################################
#Data loading and neural network image load.
################################################################################
load("~/Rotation 2/R scripts/Neural network/Neural network images/Keras.RData") #Network loading
Phosphosequences= read.csv("~\\Rotation 2\\Raw data\\Processed POIs sequences + attributes.csv")
zscore_matrix <- read.csv("~/Rotation 2/Raw data/Neural network/CDK zscore matrix.csv")
rownames(zscore_matrix)= zscore_matrix[,1]
zscore_matrix= zscore_matrix[,-1]
################################################################################
#Creating readable phosphosites
#inserts means that the phosphorylated residue has an X inserted in the +1
#position
################################################################################
Search_sequences= gsub('[[:punct:] ]+','',Phosphosequences$Raw_phosphoso)
Search_sequences= gsub('[[:digit:]]+', '', Search_sequences)
inserts=NULL
for(i in 1:nrow(Phosphosequences)){
  locations= str_locate(Phosphosequences$Full_sequence[i],Search_sequences[i])
  front= substr(Phosphosequences$Full_sequence[i],locations[1]-14,locations[1])
  back= substr(Phosphosequences$Full_sequence[i],locations[2]-14,locations[2])
  insert= paste0(front,Phosphosequences$Raw_phosphoso[i],back)
  inserts=c(inserts,insert)
  front=""
  back=""
}
inserts= gsub("0.0","",inserts)
inserts= gsub("0.9|0.8|0.7|0.6|0.5","X",inserts)
inserts= gsub('[[:punct:] ]+','',inserts)
inserts= gsub('[[:digit:]]+', '', inserts)
inserts[5]= "PMVEKQESENSCNKEEEEPVFTRQDSXNRSEKEEPVFTRQDSNRSEK" #As some sites too ambigous or contains certain (1) sites not covered by grep
inserts[6]= "KVAPPAVLNDISKKLLGPISXPPQPPSVSAWNKPISPPQPPSVSAWNK"
inserts[7]= "RQKQPRAGPIKAQKLLPDLSXPVENKIKAQKLPDLSPVENK"
inserts[13]= "RPGISTFGYNRNNKKKPYVSLAQQMAPPSXPSNSTPNSSSGSNGNDQLSKPNSSSGSNGNDQLSK"
inserts[14]= "RQKQPRAGPIKAQKLLPDLSXPVENKIKAQKLPDLSPVENK"
inserts[15]= "PSPQHTDQTEAFQKGGVPHPEDDHSQVEGPESLRPEDDHSXQVEGPESLR"
inserts[18]= "PSPQHTDQTEAFQKGGVPHPEDDHSQVEGPESLRPEDDHSXQVEGPESLR"
inserts[19]= "ERVGYREGPTVETKRRIQPQPPDEDGDHSXDKEDEQPQVVVLKHSDKEDEQPQVVVLK"
inserts[21]= "PSPSGGLSEEPAAKDDLDNRMPGLVGQEVGSXGEGPRTSSPLFNKGSGEGPRTSSPLFNK"
inserts[24]= "NPGASRDQTSPAVKQQGSXPVEPKQTSPAVKQGSPVEPK"
inserts[25]= "PGMDPAVLKAQLHKRRPEVDSXPGETPSWAPQPKVDSPGETPSWAPQPK"
inserts[33]= "PVQTRQHLNPPGGKTTSDIFGSXPVTATSRLAHPNKGSPVTATSRLAHPNK"
inserts[37]= "RPGISTFGYNRNNKKKPYVSLAQQMAPPSXPSNSTPNSSSGSNGNDQLSKPNSSSGSNGNDQLSK"
inserts[41]= "KYLCASSVGGETLDKAVCSLQKETPLPVSLPSDKTMVMEALSXLAKSSSHLSPSEEVRCTQDFLSQTQSL"
inserts[42]= "KYLCASSVGGETLDKAVCSLQKETPLPVSLPSDKTXMVMEALSLAKSSSHLSPSEEVRCTQDFLSQTQSL" #These are the extra FAM208B sites included
inserts[43]= "QDRMLCDIALWSTYGAMIPTQLPQELDFKYVMKVSXSLKKRLPEAAFRKQNYLEE"
inserts= c(inserts,"QDRMLCDIALWSTYGAMIPTQLPQELDFKYXVMKVSSLKKRLPEAAFRKQNYLEE",
           "QDRMLCDIALWSTYGAMIPTQLPQELDFKYVMKVSSXLKKRLPEAAFRKQNYLEE")
a=str_count(inserts,"X")
which(a != 1); rm(a) #Quick sanity check to ensure only one residue annotated as site
################################################################################
#Creating dataframe for prediction
################################################################################
POIs_frame=data.frame(AA1_zscore="",
                       AA2_zscore="",
                       AA3_zscore="",
                       AA4_zscore="",
                       AA5_zscore="",
                       AA6_zscore="",
                       AA7_zscore="",
                       AA8_zscore="",
                       AA9_zscore="",
                       AA10_zscore="",
                       AA11_zscore="",
                       AA12_zscore="",
                       AA13_zscore="",
                       AA14_zscore="",
                       Substrate=c(Phosphosequences$Gene_name,
                                   "FAM208B","FAM208B"))
POIs_frame$Original_phosphopeptide= c(Phosphosequences$Raw_phosphoso,"","")
POIs_frame$Cy_motif= c(Phosphosequences$Contains_Cy_CDK_motif,F,F) #as included two extra FAM208B sites.
POIs_frame$NLXXL_motif= c(Phosphosequences$Contains_NLXXL,F,F)
POIs_frame$Common_motif= c(Phosphosequences$Contains_Common,F,F)
POIs_frame$Sphase= c(Phosphosequences$Contains_S_phase,F,F)
POIs_frame$LXF= c(Phosphosequences$Contains_LXF_motif,F,F)
POIs_frame$Substrate_contains_potential_Cks_binding_site= F
POIs_frame$Substrate_contains_potential_Cks_binding_site[
  POIs_frame$Substrate %in% c("FAM207A","PRRC2C")] =T

################################################################################
#Adding amino acid symbols to data frame
################################################################################
for(i in 1:length(inserts)){
  one_up_from_target=str_locate(inserts[i],"X")
  one_up_from_target=one_up_from_target[1]
  POIs_frame$AA7_zscore[i]=substr(inserts[i],one_up_from_target-1,one_up_from_target-1)
  POIs_frame$AA6_zscore[i]=substr(inserts[i],one_up_from_target-2,one_up_from_target-2)
  POIs_frame$AA5_zscore[i]=substr(inserts[i],one_up_from_target-3,one_up_from_target-3)
  POIs_frame$AA4_zscore[i]=substr(inserts[i],one_up_from_target-4,one_up_from_target-4)
  POIs_frame$AA3_zscore[i]=substr(inserts[i],one_up_from_target-5,one_up_from_target-5)
  POIs_frame$AA2_zscore[i]=substr(inserts[i],one_up_from_target-6,one_up_from_target-6)
  POIs_frame$AA1_zscore[i]=substr(inserts[i],one_up_from_target-7,one_up_from_target-7)
  POIs_frame$AA8_zscore[i]=substr(inserts[i],one_up_from_target+1,one_up_from_target+1)
  POIs_frame$AA9_zscore[i]=substr(inserts[i],one_up_from_target+2,one_up_from_target+2)
  POIs_frame$AA10_zscore[i]=substr(inserts[i],one_up_from_target+3,one_up_from_target+3)
  POIs_frame$AA11_zscore[i]=substr(inserts[i],one_up_from_target+4,one_up_from_target+4)
  POIs_frame$AA12_zscore[i]=substr(inserts[i],one_up_from_target+5,one_up_from_target+5)
  POIs_frame$AA13_zscore[i]=substr(inserts[i],one_up_from_target+6,one_up_from_target+6)
  POIs_frame$AA14_zscore[i]=substr(inserts[i],one_up_from_target+7,one_up_from_target+7)}
letters=rownames(zscore_matrix)
################################################################################
#Ensuring empty sites covered
################################################################################
POIs_frame$AA1_zscore[POIs_frame$AA1_zscore == ""] = "X"
POIs_frame$AA2_zscore[POIs_frame$AA2_zscore == ""] = "X"
POIs_frame$AA3_zscore[POIs_frame$AA3_zscore == ""] = "X"
POIs_frame$AA4_zscore[POIs_frame$AA4_zscore == ""] = "X"
POIs_frame$AA5_zscore[POIs_frame$AA5_zscore == ""] = "X"
POIs_frame$AA6_zscore[POIs_frame$AA6_zscore == ""] = "X"
POIs_frame$AA7_zscore[POIs_frame$AA7_zscore == ""] = "X"
POIs_frame$AA8_zscore[POIs_frame$AA8_zscore == ""] = "X"
POIs_frame$AA9_zscore[POIs_frame$AA9_zscore == ""] = "X"
POIs_frame$AA10_zscore[POIs_frame$AA10_zscore == ""] = "X"
POIs_frame$AA11_zscore[POIs_frame$AA11_zscore == ""] = "X"
POIs_frame$AA12_zscore[POIs_frame$AA12_zscore == ""] = "X"
POIs_frame$AA13_zscore[POIs_frame$AA13_zscore == ""] = "X"
POIs_frame$AA14_zscore[POIs_frame$AA14_zscore == ""] = "X"
################################################################################
#Assigning Z-scores to matrix
################################################################################
letters=rownames(zscore_matrix) #Provides order of AA in matrix
fill_a_column= function(the_row){
  dict= t(zscore_matrix[,the_row]) #Defines the values a quadrant can take based on Z-scores of the corrosponding AA
  if(class(POIs_frame[,the_row]) == "numeric"){stop(... = "INVALID ROW. NOT NUMERIC")}
  for(i in 1:nrow(POIs_frame)){
    index=which(POIs_frame[i,the_row] == letters) #Gets the value assigned to the amino acid
    POIs_frame[i,the_row]= dict[index]
  }
  POIs_frame[,the_row]= as.numeric(POIs_frame[,the_row])
  POIs_frame
}
POIs_frame= fill_a_column(1)  
POIs_frame= fill_a_column(2)  
POIs_frame= fill_a_column(3)  
POIs_frame= fill_a_column(4)  
POIs_frame= fill_a_column(5)  
POIs_frame= fill_a_column(6) 
POIs_frame= fill_a_column(7)  
POIs_frame= fill_a_column(8)  
POIs_frame= fill_a_column(9)  
POIs_frame= fill_a_column(10)  
POIs_frame= fill_a_column(11)  
POIs_frame= fill_a_column(12) 
POIs_frame= fill_a_column(13)  
POIs_frame= fill_a_column(14)
POIs_frame$Cy_motif= as.numeric(POIs_frame$Cy_motif)
POIs_frame$NLXXL_motif= as.numeric(POIs_frame$NLXXL_motif)
POIs_frame$Common_motif= as.numeric(POIs_frame$Common_motif)
POIs_frame$Sphase= as.numeric(POIs_frame$Sphase)
POIs_frame$LXF= as.numeric(POIs_frame$LXF)
POIs_frame$Substrate_contains_potential_Cks_binding_site= as.numeric(POIs_frame$Substrate_contains_potential_Cks_binding_site)
write.csv(POIs_frame,"D:\\POIs_frame.csv")
