##################################################################################################
##### subtypes analysis #####
##################################################################################################
##################################################################################################
setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Subtypes")
source("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Codes/library_CIS_imagingData_12222021.R")
##################################################################################################
#subtypes 
Subtypes_analysis = subset(filterTCGA02, sampleID %in% colnames(exprFinal_mLDA))
Subtypes_analysis = as.data.frame(cbind(Subtypes_analysis$sampleID,Subtypes_analysis$PAM50mRNAnature2012))
colnames(Subtypes_analysis) = c("sampleID","PAM50mRNAnature2012")
Subtypes_analysis = Subtypes_analysis[Subtypes_analysis$PAM50mRNAnature2012!="",]
table(Subtypes_analysis$PAM50mRNAnature2012)
#Luminal A & Normal-like as 1
#Luminal B as 2
#Basal-lik as 3
#HER2-enriched as 4
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Luminal A"] = 1
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Luminal B"] = 2
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Basal-like"] = 3
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="HER2-enriched"] = 4
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Normal-like"] = 1

########DE analysis
ExprDD <- as.data.frame(exprFinal[,3:dim(exprFinal)[2]])
ExprDD <- ExprDD[,colnames(ExprDD)%in%Subtypes_analysis$sampleID,]
rownames(ExprDD) <- exprFinal$gene
temp <- as.data.frame(colnames(ExprDD))
Subtypes_analysis <- merge(temp, Subtypes_analysis, by.x ="colnames(ExprDD)", by.y = "sampleID")
#caseControlResponse <- rep(0, (dim(Surv_TCGA_10_years)[1]))
#caseControlResponse[which(Surv_TCGA_10_years$OS.time > 3650)] <- 1
library(edgeR)
group <- caseControlResponse
y <- DGEList(counts=ExprDD,group=group)
keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)

fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)
qlfT <- qlf$table
qlfT_sort <- qlfT[order(qlfT[,4]),]
qlfT_sort1<-qlfT_sort
qlfT_sort1$FC = ifelse(qlfT_sort1$logFC>0,2^(qlfT_sort1$logFC),-2^(-qlfT_sort1$logFC))
write.csv(qlfT_sort1, "survival_edgeR.csv")


#####mLDA
#subtypes 
Subtypes_analysis = subset(filterTCGA02, sampleID %in% colnames(exprFinal_mLDA))
Subtypes_analysis = as.data.frame(cbind(Subtypes_analysis$sampleID,Subtypes_analysis$PAM50mRNAnature2012))
colnames(Subtypes_analysis) = c("sampleID","PAM50mRNAnature2012")
Subtypes_analysis = Subtypes_analysis[Subtypes_analysis$PAM50mRNAnature2012!="",]
table(Subtypes_analysis$PAM50mRNAnature2012)
#Luminal A & Normal-like as 1
#Luminal B as 2
#Basal-lik as 3
#HER2-enriched as 4
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Luminal A"] = 1
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Luminal B"] = 2
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Basal-like"] = 3
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="HER2-enriched"] = 4
Subtypes_analysis$PAM50mRNAnature2012[Subtypes_analysis$PAM50mRNAnature2012=="Normal-like"] = 1

ExprDD <- as.data.frame(exprFinal[,3:dim(exprFinal)[2]])
ExprDD <- ExprDD[,colnames(ExprDD)%in%Subtypes_analysis$sampleID,]
rownames(ExprDD) <- exprFinal$gene
temp <- as.data.frame(colnames(ExprDD))
Subtypes_analysis <- merge(temp, Subtypes_analysis, by.x ="colnames(ExprDD)", by.y = "sampleID")

try2  <- mLDA.Kclass(t(ExprDD), Subtypes_analysis, t(ExprDD), Z=NULL, Z.new=NULL, K=4, tau=50, alpha=0.9, nu=150, d=2, nb=5)
saveRDS(try2,file=paste("survival_mLDA", ".rds", sep=""))
try2 <- readRDS(file=paste("survival_mLDA", ".rds", sep=""))
write.csv(as.data.frame(try2[["network.nodes"]]),"survival_mLDA.csv")
write.csv(as.data.frame(try2[["MIset"]]),"survival_mLDA_strong.csv")
survival_mLDA_network <- NULL
for (i in 1:44) {
  temp <- as.data.frame(try2[["local.network_list"]][[i]])
  temp$gene <- rownames(temp)
  temp$net <- paste("Network",i, sep = "")
  temp <- temp[,-1]
  survival_mLDA_network <- rbind(survival_mLDA_network,temp)
}
write.csv(survival_mLDA_network,"survival_mLDA_network.csv")
