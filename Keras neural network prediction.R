################################################################################
#Package install
################################################################################
if (!requireNamespace("yardstick", quietly = TRUE))
  install.packages("yardstick")
library(yardstick)
if (!requireNamespace("keras", quietly = TRUE))
  install.packages("keras")
library(keras)
################################################################################
#Loading prepared NN (Neural Network) dataframe.
################################################################################
NNdataframe= read.csv("~/Documents/Rotation_2/Raw data/Training matrix.csv")
NNdataframe$Kinase= as.numeric(NNdataframe$Kinase)
NNdataframe$Substrate_contains_potential_Cks_binding_site=
  as.numeric(NNdataframe$Substrate_contains_potential_Cks_binding_site)
NNdataframe$NLXXL_motif= as.factor(NNdataframe$NLXXL_motif)
NNdataframe$Cy_motif= as.factor(NNdataframe$Cy_motif)#
NNdataframe$Common_motif= as.factor(NNdataframe$Common_motif)
NNdataframe$Sphase= as.factor(NNdataframe$Sphase)
NNdataframe$LXF= as.factor(NNdataframe$LXF)
as.factor(data$Kinase)
columns=c("AA1_zscore","AA2_zscore","AA3_zscore","AA4_zscore","AA5_zscore",
          "AA6_zscore","AA7_zscore","AA8_zscore","AA9_zscore","AA10_zscore",
          "AA11_zscore","AA12_zscore","AA13_zscore","AA14_zscore","Kinase",
          "LXF","NLXXL_motif","Cy_motif",
          "Substrate_contains_potential_Cks_binding_site")
NNdataframe= NNdataframe[,colnames(NNdataframe) %in% columns]
motif_columns= c("LXF","Cy_motif","NLXXL_motif",
                 "Substrate_contains_potential_Cks_binding_site")
################################################################################
#Creating the training and testing set
################################################################################
set.seed(0451) #Ensure consistency between reruns
CDK1_only= which(NNdataframe$Kinase == 1)
CDK1_only_training= sample(CDK1_only,390)
`%notin%` <- Negate(`%in%`)
CDK1_only_testing= CDK1_only[CDK1_only  %notin% CDK1_only_training]
################################################################################
#Assinging the index values for the testing and training non-cdk sites
################################################################################
Otherkinases= which(NNdataframe$Kinase == 0)
Otherkinases_training=sample(Otherkinases,610)
Otherkinases_testing= which(NNdataframe$Kinase == 0)
Otherkinases_testing= Otherkinases[Otherkinases  %notin% Otherkinases_training]
Otherkinases_testing=sample(Otherkinases_testing, 610)

################################################################################
#Creating the training set
################################################################################
training_data= NNdataframe[c(CDK1_only_training,Otherkinases_training),]
testing_data= NNdataframe[c(CDK1_only_testing,Otherkinases_testing),]
################################################################################
#Defining perquisite functions
################################################################################
zscore= function(a_vector){ #Calculates standard zscore
  mean=mean(a_vector)
  sd= sd(a_vector)
  (a_vector-mean)/sd
}
standardisation= function(a_vector){
  a_vector= a_vector- min(a_vector)
  denominator= max(a_vector) - min (a_vector)
  a_vector/denominator
}
`%notin%` <- Negate(`%in%`)
hot_encode_factors= function(data_frame){
  factor_colnames= colnames(data_frame)[lapply(data_frame,class) == "factor"]
  
  for(i in 1:length(factor_colnames)){
    a_factor_column= as.factor(c(t(data_frame[factor_colnames[i]])))
    totalcategory_levels= length(levels(a_factor_column))
    levels(a_factor_column) = 1:length(a_factor_column) # one hot encode classes i.e. a matrix for the groups
    a_factor_column = to_categorical(as.integer(a_factor_column) - 1,
                                     num_classes = totalcategory_levels)
    new_col_names= NULL
    for(v in 1:totalcategory_levels){
      new_col_names= c(new_col_names,paste(factor_colnames[i],v))
    }
    
    colnames(a_factor_column)= new_col_names; rm(new_col_names)
    data_frame= cbind(data_frame,a_factor_column)
    
  }
  data_frame[,colnames(data_frame) %notin% factor_colnames]}
################################################################################
#Defining easy keras function
################################################################################
dataset=training_data
category_column_name="Kinase"
normalisation= "standardisation";
units= 100; epochs= 30
dataset= training_data; category_column_name= "Kinase"

easy_keras= function(dataset,category_column_name,normalisation= "standardisation",
                     units= 1000, epochs= 30){
  #Checks#
  if(category_column_name %notin% colnames(dataset)){
    stop("Invalid category column!")
  }
  inputs= dataset[,which(colnames(dataset) != category_column_name)]
  if(any(lapply(inputs,class) == "factor")){
    inputs= hot_encode_factors(inputs)
  }
  if(any(t(lapply(inputs,class)) != "numeric")){
    stop("Input column(s) invalid type. All inputs need to be numeric")
  }
  inputs= as.data.frame(lapply(inputs,as.numeric))
  switch(normalisation, 
         "zscore"={
           inputs= as.data.frame(lapply(inputs,zscore))},
         "standardisation"={
           inputs= as.data.frame(lapply(inputs,standardisation))},
         "none"={
           inputs= as.data.frame(lapply(inputs,as.numeric))},
         {stop("Invalid normalisation method")}) #Default outcome
  #Prepin data for next steps#
  category_column= dataset[,category_column_name]
  category_column= as.factor(category_column)
  totalcategory_column= length(levels(category_column))
  levels(category_column) = 1:length(category_column) # one hot encode classes i.e. a matrix for the groups
  category_column = to_categorical(as.integer(category_column) - 1,
                                        num_classes = totalcategory_column)
  #EDIT HERE TO MODIFY NETWORK STRUCTURE#
  model = keras_model_sequential()
  model %>% #Designing how the neural network is organised
    layer_dense(input_shape = ncol(inputs), #Dimention of the input matrix i.e. training_input size
                units = units,
                activation = "relu") %>%
    layer_dropout(rate = 0.5) %>% 
    layer_dense(units = round(ncol(inputs)/2,0), activation = 'relu') %>%
    layer_dropout(rate = 0.5) %>%
    layer_dense(units = ncol(category_column), activation = "softmax")
  model %>%
    compile( #Loss function and optmiser used to reduce errors.
      loss = "categorical_crossentropy",
      optimizer = "adagrad",
      metrics = "accuracy"
    )
  #Fitting model#
  inputs= as.matrix(inputs)
  category_column_name= as.matrix(category_column)
  fit = model %>%
    fit(
      x = inputs,
      y = category_column,
      epochs = epochs,
      batch_size = 60, 
      validation_split = 0.2
    )
  model
}


network2= easy_keras(training_data,"Kinase",epochs = 50)
network= easy_keras(training_data[,colnames(training_data) %notin% motif_columns],"Kinase",epochs = 30)

################################################################################
#Testing accuracy of neural network using the testing set
################################################################################
testing_inputs= testing_data[,colnames(training_data) != "Kinase"]
testing_categories= training_data$Kinase
testing_categories= as.factor(testing_categories)
testing_inputs = as.matrix(testing_inputs)
totalclasses= length(levels(testing_categories)) 
levels(testing_categories) = 1:length(testing_categories)
testing_categories = to_categorical(as.integer(testing_categories) - 1,
                        num_classes = totalclasses)
prediction_for_network= network %>% predict_classes(testing_inputs[,colnames(testing_inputs) %notin% motif_columns])
prediction_for_network2= network2 %>% predict_classes(testing_inputs)
################################################################################
#Defining network accuracy function
################################################################################
get_NN_accuracy= function(neural_network){
  prediction= data.frame(Prediciton= na.omit(neural_network),True_kinase= testing_data$Kinase)
  prediction$Prediciton= as.factor(prediction$Prediciton)
  prediction$True_kinase= as.factor(prediction$True_kinase)
  
  
  confusionmatrix_object= conf_mat(prediction,truth= True_kinase, 
                                   estimate = Prediciton)
  
  confusionmatrix=data.frame(
    Prediction=c("Non CDK1 site","CDK1 site",
                 "Non CDK1 site","CDK1 site"),
    'Test data'=c("Non CDK1 site","Non CDK1 site",
                  "CDK1 site","CDK1 site"),
    'Data points'=confusionmatrix_object$table[1:4],
    Outcome=c("True negative","False positive","False negative",
              "True positive"))
  confusionmatrix    
  
  confusionmatrix_object=summary(confusionmatrix_object)
  
  confusionmatrix_object= confusionmatrix_object[confusionmatrix_object$.metric %in%
                                                   c("accuracy","sens","spec","ppv","npv",
                                                     "f_meas"),]
  colnames(confusionmatrix_object)=c("Metric","Estimator","Estimate")
  confusionmatrix_object[,1]=c("Accuracy","Sensitivity","Specifiicity",
                               "Positive predictive value",
                               "Negative predictive value",
                               "F means")
  confusionmatrix_object[,-2] #Summary of models accuracy
}
################################################################################
#Obtaing neural network accuracies
################################################################################
get_NN_accuracy(prediction_for_network) #No cyclin or cks motifs
get_NN_accuracy(prediction_for_network2) #Extra parameters
################################################################################
#Prediction using intial model and confusion matrix to assess accuracy
################################################################################
POIs_frame= read.csv("~//Rotation 2//Raw data//Neural network//POI dataframe.csv")
Substrate= POIs_frame$Substrate
Sequence= POIs_frame$Original_phosphopeptide
################################################################################
#Loading prepared NN (Neural Network) dataframe.
################################################################################
POIs_frame$Substrate_contains_potential_Cks_binding_site=
  as.numeric(POIs_frame$Substrate_contains_potential_Cks_binding_site)
POIs_frame$NLXXL_motif= as.numeric(POIs_frame$NLXXL_motif)
POIs_frame$Cy_motif= as.numeric(POIs_frame$Cy_motif)#
POIs_frame$Common_motif= as.numeric(POIs_frame$Common_motif)
POIs_frame$Sphase= as.numeric(POIs_frame$Sphase)
POIs_frame$LXF= as.numeric(POIs_frame$LXF)
POIs_frame= POIs_frame[,colnames(POIs_frame) %in% columns]
POIs_frame=as.matrix(POIs_frame)
POI_prediction= network2 %>% predict_classes(POIs_frame)
POI_prediction= data.frame(Substrate= Substrate,
                       Sequence= Sequence,
                       Predition= c(POI_prediction))
POI_prediction$Predition[POI_prediction$Predition == 1]= "Cdk1 site"
POI_prediction$Predition[POI_prediction$Predition == 0]= "Not a cdk1 site"
View(POI_prediction)      
#ASSUMES HIGHEST PROBABILITY RESIDUE IS PHOSPHORYLATED
#100% CERTAIN SITES X(1) INCLUDED IN TABLE AND SEPERATE ENTRIES IN
#SAME ORDER AS THEY APPEAR IN THE SEQUENCE