################################################################################
#Package install
################################################################################
if (!requireNamespace("tidyverse", quietly = TRUE))
  install.packages("tidyverse")
library(tidyverse) 
if (!requireNamespace("recipes", quietly = TRUE))
  install.packages("recipes")
library(recipes) 
if (!requireNamespace("parsnip", quietly = TRUE))
  install.packages("parsnip")
library(parsnip) 
if (!requireNamespace("yardstick", quietly = TRUE))
  install.packages("yardstick")
library(yardstick)
if (!requireNamespace("keras", quietly = TRUE))
  install.packages("keras")
library(keras)
################################################################################
#Loading prepared NN (Neural Network) dataframe.
################################################################################
NNdataframe= read.csv("~//Rotation 2//Raw data//Neural network//NN training dataframe - Copy.csv")
NNdataframe$Kinase= as.factor(NNdataframe$Kinase)
NNdataframe$Substrate_contains_potential_Cks_binding_site=
  as.factor(NNdataframe$Substrate_contains_potential_Cks_binding_site)
NNdataframe$NLXXL_motif= as.factor(NNdataframe$NLXXL_motif)
NNdataframe$Cy_motif= as.factor(NNdataframe$Cy_motif)#
NNdataframe$Common_motif= as.factor(NNdataframe$Common_motif)
NNdataframe$Sphase= as.factor(NNdataframe$Sphase)
NNdataframe$LXF= as.factor(NNdataframe$LXF)
columns=c("AA1_zscore","AA2_zscore","AA3_zscore","AA4_zscore","AA5_zscore",
      "AA6_zscore","AA7_zscore","AA8_zscore","AA9_zscore","AA10_zscore",
      "AA11_zscore","AA12_zscore","AA13_zscore","AA14_zscore","Kinase",
      "LXF","NLXXL_motif","Cy_motif",
      "Substrate_contains_potential_Cks_binding_site")
NNdataframe= NNdataframe[,colnames(NNdataframe) %in% columns]

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
#Creating the training and testing set
################################################################################
training_data= NNdataframe[c(CDK1_only_training,Otherkinases_training),]
testing_data= NNdataframe[c(CDK1_only_testing,Otherkinases_testing),]
################################################################################
#Creating model using recipe and testing model using test data
################################################################################
the_recipe <- recipe(Kinase ~ . , data= NNdataframe) %>%
  #step_discretize(range_cols, options= list(cuts= 6)) %>% #Breaks a column into 6 ranges
  #steplog
  step_dummy(all_nominal(), -all_outcomes()) %>%
  #step_center(all_predictors(), -all_outcomes()) %>%
  #step_normalize(Age,Resting_Blood_Pressure, #Normalising nurmeric variables for modelling. ENSURE VARIABLE BEING PREDICTED ISNT INCLUDED!!!
  #               Serum_Cholesterol,
  #               Max_Heart_Rate_Achieved,
  #               ST_Depression_Exercise) %>%
prep(training_data, retain = TRUE)
train_processed_data <- juice(the_recipe)
train_processed_data <- bake(the_recipe, new_data = training_data)

################################################################################
#Creating model using recipe and testing model using test data
################################################################################
#set.seed(0451) #Ensures numbers are consistent if rerun
log_regr_hd_model <- logistic_reg(mode = "classification") %>%
  set_engine("keras") %>% 
  fit(Kinase ~ ., data = train_processed_data)
################################################################################
#Coefficients for model and their ability to predict variables
################################################################################
broom::tidy(log_regr_hd_model$fit) %>%
  arrange(desc(estimate)) %>% 
  mutate(odds_ratio = exp(estimate))

################################################################################
#Prediction using intial model and confusion matrix to assess accuracy
################################################################################
testing_data <- bake(the_recipe, new_data = testing_data)
prediction1 <- predict(log_regr_hd_model, 
                       new_data = testing_data, 
                       type     = "class")
prediction1
prediction1= cbind(na.omit(prediction1),Kinase= testing_data$Kinase)

confusionmatrix_object= conf_mat(prediction1,truth    = Kinase, 
                                 estimate = .pred_class)

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
b=confusionmatrix_object
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
  as.factor(POIs_frame$Substrate_contains_potential_Cks_binding_site)
POIs_frame$NLXXL_motif= as.factor(POIs_frame$NLXXL_motif)
POIs_frame$Cy_motif= as.factor(POIs_frame$Cy_motif)#
POIs_frame$Common_motif= as.factor(POIs_frame$Common_motif)
POIs_frame$Sphase= as.factor(POIs_frame$Sphase)
POIs_frame$LXF= as.factor(POIs_frame$LXF)
columns=c("AA1_zscore","AA2_zscore","AA3_zscore","AA4_zscore","AA5_zscore",
          "AA6_zscore","AA7_zscore","AA8_zscore","AA9_zscore","AA10_zscore",
          "AA11_zscore","AA12_zscore","AA13_zscore","AA14_zscore","Kinase",
          "Substrate_contains_potential_Cks_binding_site")
POIs_frame= POIs_frame[,colnames(POIs_frame) %in% columns]
POI_prediction= predict(log_regr_hd_model, 
        new_data = POIs_frame, 
        type     = "class")
POIs_frame= cbind(c(POI_prediction),Sequence,Substrate, POIs_frame)






















################################################################################
#Validating model using resampling
################################################################################
set.seed(0451)
hd_cv_split_objects <- NNdataframe %>%
  vfold_cv(strata = Kinase) #Resampling the data
make_cv_predictions_fcn <- function(split, id){
  analysis_tbl <- analysis(split)
  trained_analysis_recipe <- prep(the_recipe ,training = analysis_tbl)
  baked_analysis_data_tbl <- bake(trained_analysis_recipe, new_data = analysis_tbl)
  model <- logistic_reg(mode = "classification") %>%
    set_engine("glm") %>%
    fit(Kinase ~ ., data = baked_analysis_data_tbl)
  assessment_tbl <- assessment(split)
  trained_assessment_recipe <- prep(the_recipe, training = assessment_tbl)
  baked_assessment_data_tbl <- bake(trained_assessment_recipe, new_data = assessment_tbl)
  tibble("id"         = id,
         "truth"      = baked_assessment_data_tbl$Kinase,
         "prediction" = unlist(predict(model, new_data = baked_assessment_data_tbl))
  )
}
cv_predictions_tbl <- map2_df(.x = hd_cv_split_objects$splits,
                              .y = hd_cv_split_objects$id,
                              ~make_cv_predictions_fcn(split = .x, id = .y))
cv_predictions_tbl

desired_metrics <- metric_set(accuracy,
                              sens,
                              spec,
                              ppv,
                              npv,
                              f_meas)

v_metrics_long_tbl <- cv_predictions_tbl %>% 
  group_by(id) %>% 
  desired_metrics(truth = truth, estimate = prediction) 
#see results
cv_metrics_long_tbl %>% head(10) %>% knitr::kable(align = rep("c", 4))


