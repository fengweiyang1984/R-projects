
setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Datasets")

####Related to Data Import, Cleaning and Preprocessing Steps 1-10


##Load R packages:

library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)
library(survival)
library(ggplot2)
library(survminer)

data(XenaData)
write.csv(XenaData, "00_tblXenaHubInfo.csv")


####[Main Text: Step 5] Select then download target data sets from Xena Data Hubs.
##[Step 5-a] Target=RSEM expected counts provided by the UCSC toil Recompute Compendium

getOption('timeout')
options(timeout=10000)

GeneExpectedCnt_toil =XenaGenerate(subset = XenaHostNames == "toilHub") %>%
	XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
	XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")
XenaQuery(GeneExpectedCnt_toil) %>%
	XenaDownload(destdir ="./")   ###, max_try = 1L)




##[Step 5-b] Target =TCGA Clinical data.

paraCohort = "TCGA Breast Cancer"
paraDatasets ="TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"


Clin_TCGA = XenaGenerate(subset = XenaHostNames =="tcgaHub") %>%
	XenaFilter(filterCohorts = paraCohort) %>%
	XenaFilter(filterDatasets = paraDatasets)
XenaQuery(Clin_TCGA) %>%
	XenaDownload(destdir = "./")


##[Step 5-c] Target = TCGA Survival Data.

Surv_TCGA =XenaGenerate(subset = XenaHostNames == "toilHub") %>%
	XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
	XenaFilter(filterDatasets = "TCGA_survival_data")
XenaQuery(Surv_TCGA) %>%
 	XenaDownload(destdir ="./")



##[Step 5-d] Target = GTEx Phenotype data.

Pheno_GTEx = XenaGenerate(subset =XenaHostNames == "toilHub") %>%
	XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
	XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype")
XenaQuery(Pheno_GTEx) %>%
	XenaDownload(destdir ="./")




#### [Main Text: Step 6] Subset expression data to include only desired samples.

## [Step 6-a] Retrive IDs for GTEx normal samples of desired tissue type(s).
BRCA_clinicalMatrix <- read.delim("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix")
filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz")
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01))

#dim(filterGTEx01)
#[1] 19131     7


paraStudy = "GTEX" # Setting "GTEx" as the study of interest
paraPrimarySiteGTEx = "Breast" # Setting "Breast" as the primary site of interest

filterGTEx02 = subset(filterGTEx01,	study == paraStudy &  primarysite == paraPrimarySiteGTEx )

#dim(filterGTEx02)
#[1] 179   7

												
## [Step 6-b] Retrive IDs for TCGA primary tumor samples of desired histological type(s).

filterTCGA01 = fread(paraDatasets)
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01))

#dim(filterTCGA01)
#[1] 1247  194


table(filterTCGA01$primarysite)
#Breast 
#  1247 

table(filterTCGA01$histologicaltype)
#                                       Infiltrating Carcinoma NOS 
#                               6                                1 
#   Infiltrating Ductal Carcinoma   Infiltrating Lobular Carcinoma 
#                             901                              216 
#             Medullary Carcinoma            Metaplastic Carcinoma 
#                               8                                9 
#Mixed Histology (please specify)               Mucinous Carcinoma 
#                              40                               18 
#                  Other, specify 
#                              48 


paraSampleType = "Primary Tumor" #Setting "Primary Tumor" as the sample type of interest.
paraPrimarySiteTCGA = "Breast" #Setting "Breast" as the primary site of interest.


filterTCGA02 = subset(filterTCGA01,	sampletype == paraSampleType &	primarysite == paraPrimarySiteTCGA) 
#dim(filterTCGA02)
#[1] 1101  194



##[Step 6-c] Merge GTEx and TCGA sample lists. Then pull expression profiles from .gz file by IDs on merged sample list.

filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, "sample")

ExprSubsetBySamp =fread("TcgaTargetGtex_gene_expected_count.gz", select=filterExpr)

#dim(ExprSubsetBySamp)
#[1] 60498  1271



####[Main Text: Step 7] Subset expression data to include only protein coding genes

probemap =fread("zz_gencode.v23.annotation.csv", select =c(1,2))
exprALL = merge(probemap, ExprSubsetBySamp, by.x ="id", by.y = "sample")
genesPC = fread("zz_gene.protein.coding.csv")
#genesPC = read.csv("zz_gene.protein.coding.csv")
exprPC = subset(exprALL, gene %in% genesPC$Gene_Symbol)

# Remove duplicate gene symbols.
exprFinal = exprPC[!(duplicated(exprPC$gene)|
											duplicated(exprPC$gene, fromLast =TRUE)), ]

dim(exprFinal)
#[1] 18232  1272

##################################################################################################
##### case-control classification analysis
##################################################################################################
##################################################################################################
#### single gene based DE analysis
#### 
sample_vec <- sapply(colnames(exprFinal), function(x) unlist(strsplit(x, "-"))[[1]][1]) 
table(sample_vec)
#sample_vec
#gene GTEX   id TCGA 
#1  178    1 1092 
sample_vec<-sample_vec[-1:-2]
caseControlResponse <- rep(0, (dim(exprFinal)[2]-2))
caseControlResponse[which(sample_vec == "TCGA")] <- 1


ExprDD <- as.data.frame(exprFinal[,3:dim(exprFinal)[2]])
rownames(ExprDD) <- exprFinal$gene
dim(ExprDD)


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
write.csv(qlfT_sort1, "case_control_edgeR.csv")
##########################################################################################
###########################################################################################
##### detect predictive gene networks using top 200 DE genes as hub genes
#####
rownames(exprFinal) <- exprFinal[,2]
mid <- match(rownames(qlfT_sort1), rownames(exprFinal))
exprFinal_mLDA <- exprFinal[mid,]

source("library_CIS_imagingData_12222021.R")

geneName <- row.names(exprFinal_mLDA)

X <- as.matrix(t(exprFinal_mLDA))
Xsd <- apply(X, 2, sd)
length(which(Xsd==0))
#[1] 0

Xnew <- X[-1:-2,] #X[,-which(Xsd==0)]
geneName.new <- geneName[-which(Xsd==0)]
Xnew <- apply(Xnew, 2, function(x) (x-mean(x))/sd(x))


### Need to add in covariates such as age, gender,etc.

Y <- as.vector(caseControlResponse)
#system.time(try2  <- mLDA.pair(Xnew, Y, Xnew, Z=NULL, Z.new=NULL, pair=c(1,0), tau=50, alpha=0.9, nu=150, d=2, nb=5))
pred <- as.numeric(try2$PredClass)
length(which(pred!=Y))
#2

##################################################################################################
##### double-layers Cross-validation #####
##################################################################################################
##################################################################################################

#require(dplyr)
#178 GTEx 
#1092 TCGA
source("library_CIS_imagingData_12222021.R")
X=exprFinal_mLDA_cross; Y=caseControlResponse; n=ncol(exprFinal_mLDA);
sepn = 5;tau_seq = 5;nu_seq = 5
mLDA.pair.crossvalidation <-
  function(X, Y, n, sepn = sepn,tau_seq = tau_seq, nu_seq = nu_seq){
    valid_output = NULL
    temp_output = data.frame(vi=0,tau=0, nu=0, error =0)
    for (i in 1:sepn) {
      #### Preparing the training & testing dataset
      lowb1 = ceiling(178/sepn*(i-1))+1
      upb1 = min(ceiling(178/sepn*i),178)
      lowb2 = ceiling((n-178)/sepn*(i-1))+1
      upb2 = min(ceiling((n-178)/sepn*i),n-180)
      data_set1 = X[,1:178]
      data_set2 = X[,-1:-178]

      data_set1_Train = data_set1[, -lowb1:-upb1]
      data_set2_Train = data_set2[, -lowb2:-upb2]
      X_train = t(cbind(data_set1_Train,data_set2_Train))
      data_set1_Test = data_set1[, lowb1:upb1]
      data_set2_Test = data_set2[, lowb2:upb2]
      X_test = t(cbind(data_set1_Test,data_set2_Test)) 
      response_set1 = Y[1:178]
      response_set2 = Y[-1:-178]
      
      y_Train = as.vector(c(response_set1[-lowb1:-upb1],response_set2[-lowb2:-upb2]))
      y_Test =  as.vector(c(response_set1[lowb1:upb1],response_set2[lowb2:upb2]))

      #### mLDA fit ####
      for (r in 1:tau_seq) {
        tau_r = 5 * r 
        for (c in 1:nu_seq) {
          nu_c = 50 * c
          try1  <- mLDA.pair(X_train, y_Train, X_test, Z=NULL, Z.new=NULL, pair=c(1,0), tau=tau_r, alpha=0.5, nu=nu_c, d=2, nb=5)
          temp_output[1] = i
          temp_output[2] = tau_r  
          temp_output[3] = nu_c  
          temp_output[4] = length(which(try1$PredClass != y_Test))  
          valid_output = rbind(valid_output,temp_output)
          print(c(i,tau_r,nu_c))
        }
      }
    }
    return(valid_output)
  }

#5-fold double layers cross-validation
exprFinal_mLDA_cross = exprFinal_mLDA[,-1:-2]
exprFinal_mLDA_cross <- apply(exprFinal_mLDA_cross, 1, function(x) (x-mean(x))/sd(x))
exprFinal_mLDA_cross <- t(exprFinal_mLDA_cross)
resutl <-mLDA.pair.crossvalidation(X=exprFinal_mLDA_cross, Y=caseControlResponse, n=ncol(exprFinal_mLDA),
                         sepn = 5,tau_seq = 5, nu_seq = 5)


temp_output = data.frame(tau=0, nu=0, error =0);test_error_CV=NULL
for (r in 1:tau_seq) {
  tau_r = 5 * r 
  for (c in 1:nu_seq) {
    nu_c = 50 * c
    temp = subset(valid_output1,tau==tau_r&nu==nu_c)
    temp_output[1] = tau_r  
    temp_output[2] = nu_c  
    temp_output[3] = sum(temp[,4]) 
    test_error_2lCV = rbind(test_error_CV,temp_output)
  }}

################################################################
################################################################
#case-contral study using the turnning parameters decided by CV
try2  <- mLDA.pair(exprFinal_mLDA_cross, caseControlResponse, exprFinal_mLDA_cross, Z=NULL, Z.new=NULL, pair=c(1,0), tau=50, alpha=0.9, nu=150, d=2, nb=5)
saveRDS(try2,file=paste("case_contral_mLDA", ".rds", sep=""))
write.csv(as.data.frame(try2[["network.nodes"]]),"case_contral_mLDA.csv")
write.csv(as.data.frame(try2[["MIset"]]),"case_contral_mLDA_strong.csv")
#write.csv(as.data.frame(try2[["local.network_list"]][[20]]),"temp.csv")
#73 genes overlapping edgeR result
#Permutation pvalue Caculation
#setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control")
#try2 = readRDS(file=paste("case_contral_mLDA", ".rds", sep=""))
#X_All <- t(exprFinal_mLDA_cross[,colnames(exprFinal_mLDA_cross)%in%rownames(as.data.frame(try2$network.nodes))])
#Omega.r <- try2$network.Omega

################################################################
################################################################
#Permutation pvalue Caculation
setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control")
try2 = readRDS(file=paste("case_contral_mLDA", ".rds", sep=""))
All_CpG = rownames(as.data.frame(try2[["network.nodes"]]))
X_All <- exprFinal_mLDA_cross[,colnames(exprFinal_mLDA_cross)%in%All_CpG]
# scale dataset 
# for this dataset exprFinal_mLDA_cross has been normalized
#X_All = t(scale(t(X_All), center = TRUE, scale = TRUE))
All_y = as.data.frame(caseControlResponse)
X_All = cbind(All_y, X_All)
#X_All=as.data.frame(t(X_All))
Omega.r <- try2[["network.Omega"]]

ISresample_12M <- NULL
#Permutation_times = 200000
for (i in 1:200000){
  #i==1
  set.seed(i*2)
  OSD0_subsample_id <- sample(c(1:1270), 178, replace=FALSE)
  IS_M <- matrix(0, dim(X_All[,])[2], 1)
  mean12.X <- as.matrix(apply(X_All[OSD0_subsample_id,-1], 2, mean) - apply(X_All[-OSD0_subsample_id,-1], 2, mean))
  IS_M <-  Omega.r%*%mean12.X
  rownames(IS_M) <- colnames(X_All[,-1])
  ISresample_12M <- cbind(ISresample_12M, IS_M) 
  print(i)
}

X_subsample <- X_All
mean12.X <- as.matrix(apply(X_subsample[1:178,-1], 2, mean) - apply(X_subsample[-1:-178,-1], 2, mean))
IS_M <-  Omega.r%*%mean12.X

ISresample_12M = cbind(IS_M,ISresample_12M)
ISresample_12M = as.data.frame(ISresample_12M)
pvalue_CpG_sites=matrix(0, dim(X_All[,-1])[2], 2)
pvalue_CpG_sites[,1]=colnames(X_All[,-1])

ISresample_12M_A <- matrix(rep(abs(ISresample_12M[,1]), 188112), dim(ISresample_12M)[1], 188112) 
#Here abs function is making sure there is no error for further steps
ISresample_12M_B <- apply(ISresample_12M[, -1], 1, abs)
ISresample_12M_B <- apply(ISresample_12M[, -1], 2, abs)
CCC <- ISresample_12M_B-ISresample_12M_A 
ddd <- apply(CCC, 1, function(x) length(which(x>0)))/188112 
pvalue_CpG_sites[,2] = ddd
pvalue_CpG_sites = as.data.frame(pvalue_CpG_sites)
pvalue_CpG_sites$V2[pvalue_CpG_sites$V2==0] <- 10e-10
write.csv(pvalue_CpG_sites, "C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/case_control_mLDA_sig.csv")

######################################################################
###########Hotelling’s test on a selected network###################
######################################################################
nodes_name <- rownames(as.data.frame(try2[["local.network_list"]][[1]]))
XpairSel <- exprFinal_mLDA_cross[,colnames(exprFinal_mLDA_cross)%in%nodes_name]
#XpairSel = t(XpairSel)
ypair <- caseControlResponse
library(DescTools)
fit <- HotellingsT2Test(XpairSel ~ ypair)
fit$p.value
Hotelling_T2test <- Hotelling_T2test[order(Hotelling_T2test$pvalue),]

######################################################################
###################   Visualization of functional enrichment result(GSEA) ##
######################################################################
library(ggplot2)
#edgeR output
case_control_edgeR_GSEA_output_upregulated <- read_excel("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/GSEA_CC/case_control_edgeR_GSEA_output_upregulated.xlsx")
colnames(case_control_edgeR_GSEA_output_upregulated)
Upregulated = subset(case_control_edgeR_GSEA_output_upregulated, SIZE!=1)
colnames(Upregulated)[7] = "NOM_pval"
ggplot(Upregulated) +
  geom_point(aes(x = ES, y = reorder(GS, ES), color = NOM_pval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")

case_control_edgeR_GSEA_output_downregulated <- read_excel("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/GSEA_CC/case_control_edgeR_GSEA_output_downregulated.xlsx")
colnames(case_control_edgeR_GSEA_output_downregulated)
Downregulated = subset(case_control_edgeR_GSEA_output_downregulated, SIZE!=1)
colnames(Downregulated)[7] = "NOM_pval"
ggplot(Downregulated) +
  geom_point(aes(x = ES, y = reorder(GS, ES), color = NOM_pval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")

#mLDA output
case_control_mLDA_GSEA_output_upregulated <- read_excel("GSEA_CC/case_control_mLDA_GSEA_output_upregulated.xlsx")
colnames(case_control_mLDA_GSEA_output_upregulated)
Upregulated = subset(case_control_mLDA_GSEA_output_upregulated, SIZE!=1)
colnames(Upregulated)[7] = "NOM_pval"
ggplot(Upregulated) +
  geom_point(aes(x = ES, y = reorder(GS, ES), color = NOM_pval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")

case_control_mLDA_GSEA_output_downregulated <- read_excel("GSEA_CC/case_control_mLDA_GSEA_output_downregulated.xlsx")
colnames(case_control_mLDA_GSEA_output_downregulated)
Downregulated = subset(case_control_mLDA_GSEA_output_downregulated, SIZE!=1)
colnames(Downregulated)[7] = "NOM_pval"
ggplot(Downregulated) +
  geom_point(aes(x = ES, y = reorder(GS, ES), color = NOM_pval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")
