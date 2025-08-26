##################################################################################################
##### survival analysis #####
##################################################################################################
##################################################################################################

##################################################################################################
#5-years survival
Surv_TCGA_5_years = subset(TCGA_survival_data, sample %in% colnames(exprFinal_mLDA))
#patients live more than 5 years(1825 days) as events
Surv_TCGA_5_years$OS5 = 0
Surv_TCGA_5_years$OS5[Surv_TCGA_5_years$OS.time<=1825 & Surv_TCGA_5_years$OS==1] = 1
Surv_TCGA_5_years$OS.time5 = Surv_TCGA_5_years$OS.time
Surv_TCGA_5_years$OS.time5[Surv_TCGA_5_years$OS.time>1825]=1825 
#patients Phenotype to define subtype
Pheno_TCGA = subset(BRCA_clinicalMatrix, sampleID %in% Surv_TCGA_5_years$sample)
Surv_TCGA_5_years = merge(Surv_TCGA_5_years, Pheno_TCGA, by.x ="sample", by.y = "sampleID")
#colnames(Surv_TCGA_5_years)
sum(Surv_TCGA_5_years$OS5==1)/dim(Surv_TCGA_5_years)[1]

#subtypes
#AJCC_Stage 
#The staging system most often used for breast cancer is the American Joint Committee on Cancer (AJCC) TNM system.
#pathologic stage (also called the surgical stage)
km_fit <- survfit(Surv(OS.time5, OS5) ~ as.factor(AJCC_Stage_nature2012), data=Surv_TCGA_5_years)
ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Overall survival probability",
  legend = c(0.2, 0.4))
#HER2
#Breast cancer has four primary molecular subtypes, defined in large part by hormone receptors (HR)  
km_fit <- survfit(Surv(OS.time5, OS5) ~ as.factor(HER2_Final_Status_nature2012), data=Surv_TCGA_5_years)
ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Overall survival probability",legend = c(0.25, 0.2))

##################################################################################################
###### Event rates ##########
temp <- data.frame(Years=c(1:10),Events=c(20,39,69,84,99,107,120,128,133,139))
temp$rate <- round(temp$Events/1091,5)
library(ggplot2)
ggplot(data = temp, aes(x = Years, y = rate)) +
  geom_line()


####################### Take 10-years survival as thresholding to seperate 2 groups
####################### objects' survival <= 10y (3650 days) & event =1 (death) group as 0
####################### objects' survival > 10y (3650 days) group as 1
setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Survival")


########DE analysis
Surv_TCGA_10_years = subset(TCGA_survival_data, sample %in% colnames(exprFinal_mLDA))
temp = subset(Surv_TCGA_10_years,OS.time>3650)
Surv_TCGA_10_years = subset(Surv_TCGA_10_years,OS.time<=3650&OS==1)
Surv_TCGA_10_years = rbind(Surv_TCGA_10_years,temp)
ExprDD <- as.data.frame(exprFinal[,3:dim(exprFinal)[2]])
ExprDD <- ExprDD[,colnames(ExprDD)%in%Surv_TCGA_10_years$sample,]
rownames(ExprDD) <- exprFinal$gene
temp <- as.data.frame(colnames(ExprDD))
Surv_TCGA_10_years <- merge(temp, Surv_TCGA_10_years, by.x ="colnames(ExprDD)", by.y = "sample")
caseControlResponse <- rep(0, (dim(Surv_TCGA_10_years)[1]))
caseControlResponse[which(Surv_TCGA_10_years$OS.time > 3650)] <- 1
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
Surv_TCGA_10_years = subset(TCGA_survival_data, sample %in% colnames(exprFinal_mLDA))
temp = subset(Surv_TCGA_10_years,OS.time>3650)
Surv_TCGA_10_years = subset(Surv_TCGA_10_years,OS.time<=3650&OS==1)
Surv_TCGA_10_years = rbind(Surv_TCGA_10_years,temp)
exprFinal_mLDA_10ysurvival <- exprFinal_mLDA_cross[rownames(exprFinal_mLDA_cross)%in%Surv_TCGA_10_years$sample,]
temp <- as.data.frame(rownames(exprFinal_mLDA_10ysurvival))
Surv_TCGA_10_years <- merge(temp, Surv_TCGA_10_years, by.x ="rownames(exprFinal_mLDA_10ysurvival)", by.y = "sample")
caseControlResponse <- rep(0, (dim(Surv_TCGA_10_years)[1]))
caseControlResponse[which(Surv_TCGA_10_years$OS.time > 3650)] <- 1

try2  <- mLDA.pair(exprFinal_mLDA_10ysurvival, caseControlResponse, exprFinal_mLDA_10ysurvival, Z=NULL, Z.new=NULL, pair=c(1,0), tau=50, alpha=0.9, nu=150, d=2, nb=5)
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
