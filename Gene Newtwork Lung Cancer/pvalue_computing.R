setwd("C:\\Users\\fengw\\OneDrive\\??????\\Paper work\\Results\\v2.0")
################################################################
################################################################
#Select the combined component(biomarker)
marker_all = NULL
for (cc in 1:22) {
  #cc=1
  chr_n=paste("chr", cc, sep="")
  try1 = readRDS(file=paste("MI_CpG_chr", cc, ".rds", sep=""))
  marker_chr = illuminaMethyl450_hg38_GDC[illuminaMethyl450_hg38_GDC$Chromsome==chr_n,]
  methylation_data_10k = methylation_data1[rownames(methylation_data1)%in%marker_chr$Marker,]
  All_CpG = rownames(methylation_data_10k[try1$screenset,])
  temp = as.data.frame(All_CpG)
  marker_all = rbind(marker_all,temp)
}

#for (cc in 1:22) {
  cc=21
  chr_n=paste("chr", cc, sep="")
  try1 = readRDS(file=paste("MI_CpG_chr", cc, ".rds", sep=""))
  temp = as.data.frame(rownames(as.data.frame(try1[["local.network_list"]][[3]])))
  #marker_all = rbind(marker_all,temp)
#}
write.csv(temp,"temp.csv")
################################################################
################################################################
#Permutation pvalue Caculation
Z1=as.data.frame(t(Z))
Z1=as.data.frame(t(Z1[with(Z1, order(OSD)),]))
All_y <- Z1[1,]
X_All <- methylation_data1[rownames(methylation_data1)%in%marker_all$All_CpG,]
X_All[which(X_All==0, arr.ind=TRUE)] <- 0.01
X_All[is.na(X_All)]=0.01
X_All[which(X_All=="na", arr.ind=TRUE)]=0.01
X_All[which(X_All=="N/A", arr.ind=TRUE)]=0.01
X_All <- log(X_All)
# scale dataset 
X_All = t(scale(t(X_All), center = TRUE, scale = TRUE))
X_All = rbind(All_y, X_All)
X_All=as.data.frame(t(X_All))
X_All=as.data.frame(t(X_All[with(X_All, order(OSD)),]))
Z0 = Z1[-1,]
Omega.r <- my.inv(cov(t(X_All[-1,])))
Omega.Z <- qr.solve(cov(t(Z1[-1,])))

ISresample_12M <- NULL
for (i in 1:200000){
  #i==1
  set.seed(i*2)
  OSD0_subsample_id <- sample(c(1:492), 246, replace=FALSE)
  #OSD1_subsample_id 
  #X_subsample <- as.data.frame(t(X_All[,c(OSD0_subsample_id,OSD1_subsample_id)]))
  #Z_subsample <- Z1[,c(OSD0_subsample_id,OSD1_subsample_id)]
  temp=Z1[1,]
  Z_subsample <- apply(Z1[-1,], 1, function(x) (x-mean(x))/sd(x))
  Z_subsample = as.data.frame(t(Z_subsample))
  Z1 = as.data.frame(t(rbind(temp,Z_subsample)))
  IS_M <- matrix(0, dim(X_All[-1,])[1], 1)
  mean12.X <- as.matrix(apply(X_All[-1,OSD0_subsample_id], 1, mean) - apply(X_All[-1,-OSD0_subsample_id], 1, mean))
  IS_M <-  Omega.r%*%mean12.X
  row.names(IS_M) <- rownames(X_All[-1,])
  ISresample_12M <- cbind(ISresample_12M, IS_M)  
}

X_All_1 = t(X_All)
X_subsample <- as.data.frame(t(X_All))
mean12.X <- as.matrix(apply(t(X_subsample[X_subsample$OSD==0,-1]), 1, mean) - apply(t(X_subsample[which(X_subsample$OSD==1),-1]), 1, mean))
IS_M <-  Omega.r%*%mean12.X
Z_subsample <- as.data.frame(t(Z))
temp <- as.matrix((apply(t(Z_subsample[Z_subsample$OSD==0,-1]),1,mean) - apply(t(Z_subsample[which(Z_subsample$OSD==1),-1]),1,mean))%*%Omega.Z)
#IS_M <- abs(IS_M)+ abs(temp[1])*1+abs(temp[2])*2
ISresample_12M = cbind(IS_M,ISresample_12M)
ISresample_12M = as.data.frame(ISresample_12M)
pvalue_CpG_sites=matrix(0, dim(X_All[,])[1], 2   )
pvalue_CpG_sites[,1]=rownames(X_All)
pvalue_CpG_sites <- pvalue_CpG_sites[-1,]
#Matrix computing save a lot of time
ISresample_12M_A <- matrix(rep(abs(ISresample_12M[,1]), 200000), dim(ISresample_12M)[1], 200000) 
ISresample_12M_B <- abs(ISresample_12M[, -1]) 
ISresample_12M_B <- apply(ISresample_12M[, -1], 1, abs)
ISresample_12M_B <- apply(ISresample_12M[, -1], 2, abs)
CCC <- ISresample_12M_B-ISresample_12M_A 
ddd <- apply(CCC, 1, function(x) length(which(x>0)))/200000 
pvalue_CpG_sites[,2] = ddd

dim(ISresample_12M)
ISresample_12M <- ISresample_12M[,-1]
rm(ISresample_12M_A)
rm(ISresample_12M_B)
rm(CCC)
rm(ddd)

######################################################################
###########Multiple test by benjamini-hochberg procedure##############
######################################################################
pvalue_CpG_sites <- as.data.frame(pvalue_CpG_sites)
pvalue_CpG_sites <- pvalue_CpG_sites[order(pvalue_CpG_sites$V2),]
colnames(pvalue_CpG_sites) <- c("CpGs", "pval_raw")
#require("sgof")
#temp1 = BH(temp)
pvalue_CpG_sites$pval_adj=p.adjust(pvalue_CpG_sites$pval_raw, method = c("BH"), n = nrow(pvalue_CpG_sites))
write.csv(pvalue_CpG_sites, file="permutation_pvalue_CpG_sites.csv")

######################################################################
###########Hotelling’s test on a selected network###################
######################################################################

#try1$MIset
#network_Omega2 <- network_Omega1[row.names(network_Omega1)=="cg22301261",] 


Adj <- abs(Omega.r)
colnames(Adj) <- node_r
rownames(Adj) <- node_r

vis=rep(0, dim(Adj)[1])
#source("library_CIS_imagingData_08252020.R")
tmp <- con_node(id[1], Adj, visited=vis, m=10, depth=0,  connect=NULL, path=NULL)$connect

node_r_new <- rownames(Adj)[tmp]
Adj_new <- Adj[tmp, tmp]
length(node_r_new)   #### network size

tmp0 <- match(Sele_node, node_r_new)
tmp0 <- tmp0[which(!is.na(tmp0))]
Sele_node_new <- node_r_new[tmp0]   
Sele_node_new         #### Selected voxels in the network
length(Sele_node_new) #### Number of selected voxels in the network

Sele_marg_id <- match(Sele_node_new, rownames(log_methylation_data1))
tmp1 <- match(Sele_node_new, node_r)
Omega.r_sele <- Omega.r[tmp1, tmp1]

cc=22
chr_n=paste("chr", cc, sep="")
try1 = readRDS(file=paste("MI_CpG_chr", cc, ".rds", sep=""))
network_Omega <- try1$network.Omega
network_Omega1 <- as.data.frame(network_Omega[,row.names(network_Omega)=="cg14815005"])
network_Omega1 <- subset(network_Omega1, network_Omega1[,1]!=0)
local_network_id <- rownames(network_Omega1)
methylation_data_10k = methylation_data1[rownames(methylation_data1)%in%local_network_id,]
methylation_data_10k[which(methylation_data_10k==0, arr.ind=TRUE)] <- 0.01
# Correct the cell whose value is NA/na/N/a
methylation_data_10k[is.na(methylation_data_10k)]=0.01
methylation_data_10k[which(methylation_data_10k=="na", arr.ind=TRUE)]=0.01
methylation_data_10k[which(methylation_data_10k=="N/A", arr.ind=TRUE)]=0.01
log_methylation_data_10k <- log(methylation_data_10k)
log_methylation_data_10k = t(scale(t(log_methylation_data_10k), center = TRUE, scale = TRUE))
log_methylation_data_10k = rbind(abc,log_methylation_data_10k)

#library(DescTools)
XpairSel <- log_methylation_data_10k[-1,]
XpairSel = t(XpairSel)
ypair <- t(log_methylation_data_10k[1,])

HotellingsT2Test(XpairSel ~ ypair)

Hotelling_T2test <- Hotelling_T2test[order(Hotelling_T2test$pvalue),]
Hotelling_T2test$pval_adj=p.adjust(Hotelling_T2test$pvalue, method = c("BH"), n = 43)
write.csv(Hotelling_T2test, file="Hotelling_T2test.csv")
######################################################################
###################    survival kaplan meier      ####################
######################################################################

library(survival)
library(ggplot2)
library(survminer)
par(mfrow = c(2, 2))
#overall survival
km_fit <- survfit(Surv(OS.time, OS) ~ OSD_mLDA, data=TCGA.LUAD.survival,conf.int=.95, 
                  type=c("kaplan-meier"))
surv_pvalue(km_fit)
ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Overall survival probability", 
  pval = TRUE,
  #pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE)
#PFI, Progression Free Interval
km_fit <- survfit(Surv(OS.time) ~ OSD_mLDA, data=TCGA.LUAD.survival,conf.int=.95, 
                  type=c("kaplan-meier"))
surv_pvalue(km_fit)
ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Progression Free Interval", 
  #pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE)
#DFI, Disease Free Interval
km_fit <- survfit(Surv(OS.time,OS==1) ~ OSD_mLDA, data=TCGA.LUAD.survival,conf.int=.95, 
                  type=c("kaplan-meier"))
surv_pvalue(km_fit)
ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Disease Free Interval", 
  #pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE)
#DSS, Disease Specific Survival
km_fit <- survfit(Surv(OS.time,OS==0) ~ OSD_mLDA, data=TCGA.LUAD.survival,conf.int=.95, 
                  type=c("kaplan-meier"))
surv_pvalue(km_fit)
ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Disease Specific Survival", 
  #pval = TRUE,             # show p-value of log-rank test.
  conf.int = TRUE)

######################################################################
###################   Ttest for CpGs      ####################
######################################################################
Ttest_result <- data.frame(Tvalue=rep(0,  dim(methylation_data1)[1]))

methylation_data1_ttest <- rbind(Z[1,],methylation_data1)
methylation_data1_ttest <- as.data.frame(t(methylation_data1_ttest))
Ttest_result$meandiff <- abs(apply(t(methylation_data1_ttest[methylation_data1_ttest$OSD==0,-1]), 1, mean) - 
                        apply(t(methylation_data1_ttest[methylation_data1_ttest$OSD==1,-1]), 1, mean))

Ttest_result$deviation<- (apply(t(methylation_data1_ttest[methylation_data1_ttest$OSD==0,-1]), 1, sd)^2 + 
                               apply(t(methylation_data1_ttest[methylation_data1_ttest$OSD==1,-1]), 1, sd)^2)*245/490

Ttest_result$Tvalue <- Ttest_result$meandiff/sqrt(Ttest_result$deviation*2/246)
Ttest_result$pvalue <- pt(Ttest_result$Tvalue, df=490, lower.tail = F)*2
rownames(Ttest_result) <- rownames(methylation_data1)
Ttest_result$logrank <- -log2(Ttest_result$pvalue)
Ttest_result <- Ttest_result[order(Ttest_result[,5]),]
#Ttest_result$gene <- Annotation_CpGs[Annotation_CpGs$gene%in%rownames(Ttest_result),2]
write.csv(Ttest_result,"Ttest_result.csv")
write.csv(Annotation_CpGs,"Annotation_CpGs.csv")

######################################################################
###################   Visualization of functional enrichment result ##
######################################################################
#library(ggplot2)
ggplot(Upregulated) +
  geom_point(aes(x = ES, y = reorder(GS, ES), color = NOM.p.val, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")

