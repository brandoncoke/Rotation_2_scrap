################################################################################
#Package install
################################################################################
if("BiocManager" %in% rownames(installed.packages()) == FALSE){
  install.packages("BiocManager")}
library(BiocManager)
if("depmap" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("depmap")}
library(depmap)

################################################################################
#Package install
################################################################################
















# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("depmap")