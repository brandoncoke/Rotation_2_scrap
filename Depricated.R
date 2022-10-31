################################################################################
#Package install
################################################################################
if("neuralnet" %in% rownames(installed.packages()) == FALSE){
  install.packages("neuralnet")}
library(neuralnet)
################################################################################
#Package install
################################################################################
NNdataframe= read.csv("D:\\Training matrix.csv")
################################################################################
#Creating the training and testing set
################################################################################
set.seed(0451) #Ensure consistency between reruns
CDK1_only= which(NNdataframe$Kinase == 1)
Otherkinases= which(NNdataframe$Kinase == 0)
Otherkinases=sample(Otherkinases, 780)
training_data= NNdataframe[c(CDK1_only,Otherkinases),]
  
  
  
################################################################################
#Creating network
################################################################################
Network= neuralnet(Kinase~AA1_zscore+AA2_zscore+AA3_zscore+AA4_zscore+AA5_zscore+
                     AA6_zscore+AA7_zscore+AA8_zscore+AA9_zscore+AA10_zscore+
                     AA11_zscore+AA12_zscore+AA13_zscore+AA14_zscore+Cy_motif+
                     NLXXL_motif+Common_motif+Sphase+LXF+Substrate_contains_potential_Cks_binding_site,
                   data=training_data,
                   linear.output = FALSE,
                   hidden=c(7,5,3),
                   stepmax = 5e+5)
Test_prediction <- predict(Network, testing_data)#ENSURE training set has same variable names!!!
Test_prediction[Test_prediction > 0.5] = 1
Test_prediction[Test_prediction < 0.5] =0
#Predictions are floats 0 to 1. if >0.5 then 1
results= table(testing_data$Kinase, Test_prediction) #get an idea of accuracy in redumentary confusion matrix
results
plot(Network)
Sensitivity= results[2,2]/sum(results[,2]); Sensitivity
Specificity= results[1,1]/sum(results[,1]); Specificity
save.image(file= "Network.RData")
