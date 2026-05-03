# Task: Legacy GEO helper script retained from the original project; not part of the main PLZ-HF integration workflow.
# Manuscript mapping: archival helper only.

#rm(list=ls())
#gc()
options(stringsAsFactors = F)  
# suppressMessages(library(GEOquery))
library(GEOquery)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(limma)
library(sva)


## download data
gset = getGEO('', destdir=".",getGPL = F)
gset = getGEO('GSE116959', destdir=".", AnnotGPL = F, getGPL = F)  

exp<-exprs(gset[[1]]) 

pdata<-pData(gset[[1]]) 
## obtain group info
pdata = pData(gset)
