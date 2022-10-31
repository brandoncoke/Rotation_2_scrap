################################################################################
#Package install
################################################################################
if("readxl" %in% rownames(installed.packages()) == FALSE){  
  install.packages("readxl")}
library(readxl)
the_plate= read_xlsx("D:\\Sheet1.xlsx",
                     range= anchored("B14", dim = c(8, 12), col_names = T))
################################################################################
#Importing the 96 well into R. First file selection
#Second IF YOU KNOWN THE FILE LOCATION
################################################################################
blanks=as.numeric(the_plate[1,10:12])
blank_average= mean(blanks)
concentrations=c(125,250,500,750,1000,1500,2000) #Concentration is in mg ml
the_plate= the_plate-blank_average #This is the normalisation step
standard_dataframe= data.frame(row.names= concentrations,
                               replicate_1= the_plate[2:8,10],
                               replicate_2= the_plate[2:8,11],
                               replicate_3= the_plate[2:8,12])
colnames(standard_dataframe)= c("replicate_1","replicate_2",
                                "replicate_3")


mean_absorbance= apply(standard_dataframe[1:3], 1, mean, na.rm = TRUE)
sd_absorbance= apply(standard_dataframe[1:3], 1, sd, na.rm= T)
standard_dataframe= cbind(standard_dataframe,sd_absorbance)
################################################################################
#Importing the 96 well into R. First file selection
#Second IF YOU KNOWN THE FILE LOCATION
################################################################################
plot(concentrations,mean_absorbance)
#arrows(x0=concentrations, y0=mean_absorbance-standard_error, #Error bars.
#       x1=concentrations, y1=mean_absorbance+standard_error, code=3, angle=90)
Standard <- lm(mean_absorbance ~ concentrations)
abline(Standard)
equation= summary(Standard)
c=equation$coefficients[1]
gradient= equation$coefficients[2]
################################################################################
#Obtaining protein concentration
################################################################################
get_my_concentration= function(numeric_list,care_about_outliers= T){
  if(sd(numeric_list,na.rm = T)>0.5 && care_about_outliers || #I know two booleans bad practice 
     class(numeric_list) != "numeric"){
    print("TOO MUCH VARIATION IN SAMPLE");
    stop()}
  average_absorbance= mean(numeric_list,na.rm = T)
  sample_conc= (average_absorbance - c) / gradient #In ug. AFTER 10/3 DILUTION.
  sample_conc= sample_conc * (10/3) #Final conc in supernatant ug ml-1
  sample_conc= sample_conc / 1000 #Getting sample conc in ug ul-1
  if(sample_conc<0){
    print("invalid concentration")
  }
  else{sample_conc}
}
sample1= get_my_concentration(as.numeric(the_plate[8,7:9])) #First sample
sample2= get_my_concentration(as.numeric(the_plate[7,7:9])) #Second sample
sample3= get_my_concentration(as.numeric(the_plate[6,7:9])) #Third sample
sample4= get_my_concentration(as.numeric(the_plate[5,7:9])) #Forth sample
paste("Sample 1's concentration is", sample1, "sample 2's concentration is",
      sample2, "sample 3's concentration is", sample3,
      "sample 4's concentration", sample4)
################################################################################
#Volume needed for 25 ug per well after SDS dilution
################################################################################
50 / sample1
50 / sample2
50 / sample3
50 / sample4
################################################################################
#Warning messages
################################################################################
if(any(sd_absorbance > 0.3)){
if(sd(blanks)>0.3){print("TOO MUCH VARIATION IN THE BLANKS!!!")}
  print("YUGE VARIANTION IN ONE OR MORE OF THE POINTS"); View(standard_dataframe)
}
