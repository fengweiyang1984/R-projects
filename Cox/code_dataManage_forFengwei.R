
setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Dataset")

####Related to Data Import, Cleaning and Preprocessing Steps 1-10


##Load R packages:

library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)
library(survival)
library(ggplot2)
library(survminer)

####Load normalized gene expression data
ExprSubsetBySamp =fread("TCGA.OV.sampleMap_HiSeqV2_PANCAN.gz")
exprALL <- as.data.frame(t(ExprSubsetBySamp))
colnames(exprALL)<- exprALL[1,]; exprALL<-exprALL[-1,]

####Load survival data
survival_OV_survival <- read.delim("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Dataset/survival_OV_survival.txt")
survival_OV <- subset(survival_OV_survival, sample %in% rownames(exprALL))
survival_OV <- survival_OV[,1:4]
survival_OV[is.na(survival_OV)]<-0
