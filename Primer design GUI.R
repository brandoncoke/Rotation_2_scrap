if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)
if("openPrimeR" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("openPrimeR")}
library(openPrimeR)
if("openPrimeRui" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("openPrimeRui")}
library(openPrimeRui)
if("ViennaRNA" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("ViennaRNA")}
library(ViennaRNA)
if("OligoArrayAux" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("OligoArrayAux")}
library(OligoArrayAux)
if("MAFFT" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("MAFFT")}
library(MAFFT)
if("MELTING" %in% rownames(installed.packages()) == FALSE){
  BiocManager::install("MELTING")}
library(MELTING)





library(openPrimeR)
if (interactive()) {
  startApp()
}

BiocManager::install("ViennaRNA")