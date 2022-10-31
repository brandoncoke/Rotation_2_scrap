################################################################################
#Package install
################################################################################
if("neuralnet" %in% rownames(installed.packages()) == FALSE){
  install.packages("neuralnet")}
library(neuralnet)
################################################################################
#Normalisation functions
################################################################################
zscore= function(x,a_vector){ #Calculates standard zscore
  mean=mean(a_vector)
  sd= sd(a_vector)
  (x-mean)/sd
}
scalingtorange= function(x,a_vector){
  (x-max(a_vector))/(max(a_vector)-min(a_vector))
}
feature_clipping= function(x,a_vector){ #If output= F then outlier. Base do on zscore.
  x_zscore= zscore(x,a_vector)
  if(x_zscore > 2 || x_zscore < -2){
    output=FALSE
  } else{
    output=TRUE}
  output
}
#Alternatively log scale
################################################################################
#Read csv and simplication of columns
################################################################################
NNdataframe= read.csv("~//Rotation 2//Raw data//Neural network//CDK z-scores for training extra columns.csv")
CDK_sites= which(grepl("CDK",NNdataframe$Kinase))
CDKkinases= NNdataframe[CDK_sites,]
Otherkinases= NNdataframe[-CDK_sites,]
CDKkinases$Kinase= "CDK"
Otherkinases$Kinase= "Other"
NNdataframe= rbind(CDKkinases,Otherkinases)
NNdataframe$Kinase[NNdataframe$Kinase == "CDK"] = 1
NNdataframe$Kinase[NNdataframe$Kinase == "Other"] = 0
NNdataframe$AA1_zscore= as.numeric(NNdataframe$AA1_zscore)
NNdataframe$AA2_zscore= as.numeric(NNdataframe$AA2_zscore)
NNdataframe$AA3_zscore= as.numeric(NNdataframe$AA3_zscore)
NNdataframe$AA4_zscore= as.numeric(NNdataframe$AA4_zscore)
NNdataframe$AA5_zscore= as.numeric(NNdataframe$AA5_zscore)
NNdataframe$AA6_zscore= as.numeric(NNdataframe$AA6_zscore)
NNdataframe$AA7_zscore= as.numeric(NNdataframe$AA7_zscore)
NNdataframe$AA8_zscore= as.numeric(NNdataframe$AA8_zscore)
NNdataframe$AA9_zscore= as.numeric(NNdataframe$AA9_zscore)
NNdataframe$AA10_zscore= as.numeric(NNdataframe$AA10_zscore)
NNdataframe$AA11_zscore= as.numeric(NNdataframe$AA11_zscore)
NNdataframe$AA12_zscore= as.numeric(NNdataframe$AA12_zscore)
NNdataframe$AA13_zscore= as.numeric(NNdataframe$AA13_zscore)
NNdataframe$AA14_zscore= as.numeric(NNdataframe$AA14_zscore)
NNdataframe$Kinase= as.numeric(NNdataframe$Kinase)



################################################################################
#Creating the training and testing set
################################################################################
set.seed(0451) #Ensure consistency between reruns
trainset= sample(nrow(NNdataframe), .8*nrow(NNdataframe))
training_data= NNdataframe[trainset,]
testing_data= NNdataframe[-trainset, ] #Basically removes data points not in training set
rm(trainset)

################################################################################
#Creating network
################################################################################
Network= neuralnet(Kinase~AA1_zscore+AA2_zscore+AA3_zscore+AA4_zscore+Cy_motif+
                     NLXXL_motif+G2_motif+AA5_zscore+AA6_zscore+AA7_zscore+
                     AA8_zscore+AA9_zscore+AA10_zscore+AA11_zscore+AA12_zscore+
                     AA13_zscore+AA14_zscore,
                   data=training_data,
                   linear.output = FALSE,
                   hidden=c(5,3),
                   stepmax = 1e+6)
Test_prediction <- predict(Network, testing_data)#ENSURE training set has same variable names!!!
Test_prediction[Test_prediction > 0.5] = 1
Test_prediction[Test_prediction < 0.5] =0
#Predictions are floats 0 to 1. if >0.5 then 1
results= table(testing_data$Kinase, Test_prediction) #get an idea of accuracy in redumentary confusion matrix
results
plot(Network)
Sensitivity= results[2,2]/sum(results[,2]); Sensitivity
Specificity= results[1,1]/sum(results[,1]); Specificity
  
Test_prediction <- predict(Network, NNdataframe)#ENSURE training set has same variable names!!!
Test_prediction[Test_prediction > 0.5] = 1
Test_prediction[Test_prediction < 0.5] =0
POIs= read.csv("~//Rotation 2//Raw data//Neural network//POI z-scores with extra columns.csv")
Prediction= predict(Network, POIs)
Probable_cdk_site=Prediction
Probable_cdk_site[Prediction > 0.5] = 1
Probable_cdk_site[Prediction < 0.5] = 0

Predictions= data.frame(Substrate= POIs$Substrate,
                        Phosphosite=POIs$Original_phosphopeptide,
                        Raw_prediction=Prediction,
                        Probable_cdk_site=Probable_cdk_site)
View(Predictions)
  
  
  


