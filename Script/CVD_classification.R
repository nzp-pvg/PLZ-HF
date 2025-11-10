#rm(list=ls())
#gc()
#避免自动将字符串转换为R语言因子
options(stringsAsFactors = F)  
##加载R包
# suppressMessages(library(GEOquery))
library(GEOquery)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tidyr)
library(limma)
library(sva)

##set work path
setwd("~/Science/CVD/CVD_MS_1/data/GEO")

## download data
gset = getGEO('', destdir=".",getGPL = F)
#设置后两个为T比较耗时，而且作用不大
gset = getGEO('GSE116959', destdir=".", AnnotGPL = F, getGPL = F)  

#exp即为表达矩阵
exp<-exprs(gset[[1]]) 

## 获取ExpressionSet对象，包括的表达矩阵和分组信息gset=gset[[1]]
#使用pData()函数提取临床信息
pdata<-pData(gset[[1]]) 
## obtain group info
pdata = pData(gset)