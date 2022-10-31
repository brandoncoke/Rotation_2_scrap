################################################################################
#Package install
################################################################################

################################################################################
#Reading csvs
################################################################################
location_of_correlation_csv= "~\\Rotation 2\\Raw data\\Correlations"
files=list.files(location_of_correlation_csv)
FAM207A= read.csv(paste0(location_of_correlation_csv,"\\",files[1]))
FAM207A= FAM207A[,-c(1,11)] #Getting rid of pointless index columns
colnames(FAM207A)[c(1,4,7)]= "Gene" #This enables merge to work
in_all_three_merge= function(dataframe){
  first_data_frame=dataframe[,1:3]
  second_data_frame=dataframe[,4:6]
  third_data_frame=dataframe[,7:9]
output= merge(first_data_frame,second_data_frame)
merge(output,third_data_frame)}
FAM207A_merged= in_all_three_merge(FAM207A)
write.csv(FAM207A_merged,"~\\FAM207A.csv")
################################################################################
#Reading JPT2
################################################################################
JPT2= read.csv(paste0(location_of_correlation_csv,"\\",files[2]))
JPT2= JPT2[,-c(1,11)] #Getting rid of pointless index columns
colnames(JPT2)[c(1,4,7)]= "Gene" #This enables merge to work
JPT2_merged= in_all_three_merge(JPT2)
write.csv(JPT2_merged,"~\\JPT2.csv")
################################################################################
#Reading KIAA1143
################################################################################
KIAA1143= read.csv(paste0(location_of_correlation_csv,"\\",files[3]))
KIAA1143= KIAA1143[,-c(1,11)] #Getting rid of pointless index columns
colnames(KIAA1143)[c(1,4,7)]= "Gene" #This enables merge to work
KIAA1143_merged= in_all_three_merge(KIAA1143)
write.csv(KIAA1143_merged,"~\\KIAA1143.csv")
################################################################################
#Reading KIAA1671 
################################################################################
KIAA1671 = read.csv(paste0(location_of_correlation_csv,"\\",files[4]))
KIAA1671 = KIAA1671[,-c(1,11)] #Getting rid of pointless index columns
colnames(KIAA1671)[c(1,4,7)]= "Gene" #This enables merge to work
KIAA1671_merged= in_all_three_merge(KIAA1671)
write.csv(KIAA1671_merged,"~\\KIAA1671.csv")
################################################################################
#Reading PDXDC1 
################################################################################
PDXDC1= read.csv(paste0(location_of_correlation_csv,"\\",files[5]))
PDXDC1= PDXDC1[,-c(1,11)] #Getting rid of pointless index columns
colnames(PDXDC1)[c(1,4,7)]= "Gene" #This enables merge to work
PDXDC1_merged= in_all_three_merge(PDXDC1)
write.csv(PDXDC1_merged,"~\\PDXDC1.csv")
################################################################################
#Reading PRRC2C
################################################################################
PRRC2C= read.csv(paste0(location_of_correlation_csv,"\\",files[6]))
PRRC2C= PRRC2C[,-c(1,11)] #Getting rid of pointless index columns
colnames(PRRC2C)[c(1,4,7)]= "Gene" #This enables merge to work
PRRC2C_merged= in_all_three_merge(PRRC2C)
write.csv(PRRC2C_merged,"~\\PRRC2C.csv")
################################################################################
#Reading RBMS2
################################################################################
RBMS2= read.csv(paste0(location_of_correlation_csv,"\\",files[7]))
RBMS2= RBMS2[,-c(1,11)] #Getting rid of pointless index columns
colnames(RBMS2)[c(1,4,7)]= "Gene" #This enables merge to work
RBMS2_merged= in_all_three_merge(RBMS2)
write.csv(RBMS2_merged,"~\\RBMS2.csv")
