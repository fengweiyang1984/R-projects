setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control")
setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Other paper/Breast Cancer/Results/Case_control")
source("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Codes/library_CIS_imagingData_122022.R")
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

library( "DESeq2" )
metaData <- as.data.frame(colnames(ExprDD))
metaData$TC <- caseControlResponse
metaData$TC[metaData$TC==0] <- "control"
metaData$TC[metaData$TC==1] <- "case"
ExprDD1 <- cbind(rownames(ExprDD),ExprDD)
dds <- DESeqDataSetFromMatrix(countData=ExprDD1, colData=metaData, design=~TC, tidy = TRUE)
dds <- DESeq(dds)
resMFType <- results(dds)
#resLRT_Sig = subset(resLRT, pvalue<0.05)
resLRTOrdered <- resMFType[order(resMFType$pvalue),]
write.csv(as.data.frame(resLRTOrdered), file="case_control_Deseq.csv")
##########################################################################################
###########################################################################################
##### detect predictive gene networks using top 200 DE genes as hub genes
#####
rownames(exprFinal) <- exprFinal[,2]
mid <- match(rownames(qlfT_sort1), rownames(exprFinal))
exprFinal_mLDA <- exprFinal[mid,]


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
#case-contral study using the mLDA by the strong signals from edgeR
library(edgeR)
library( "DESeq2" )
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
qlfT_sort$fdr <- p.adjust(qlfT_sort[,4], method ="fdr")
edgeR.genes <- rownames(qlfT_sort)
length(edgeR.genes)
write.csv(as.data.frame(qlfT_sort), file="case_control_edgeR.csv")

metaData <- as.data.frame(colnames(ExprDD))
metaData$TC <- caseControlResponse
ExprDD1 <- cbind(rownames(ExprDD),ExprDD)
dds <- DESeqDataSetFromMatrix(countData=ExprDD1, colData=metaData, design=~TC, tidy = TRUE)
dds <- DESeq(dds)
resMFType <- results(dds)
#summary(resLRT) 
#resLRT_Sig = subset(resLRT, pvalue<0.05)
resLRTOrdered <- resMFType[order(resMFType$pvalue),]
write.csv(as.data.frame(resLRTOrdered), file="case_control_Deseq.csv")

#edgeR.genes <- sapply(rownames(qlfT_sort), function(x) unlist(strsplit(x, "\\."))[[1]][1])
edgeR.genes <- rownames(qlfT_sort)
edgeR.set <- match(edgeR.genes, colnames(exprFinal_mLDA_cross))
which(is.na(edgeR.set))
Xnew <- exprFinal_mLDA_cross[,edgeR.set]
#### Use the 50 most significant genes from edgeR as the strong gene set. The input is the gene positions (column numbers) in the input Xnew design matrix.  
match(c("BRCA1","BRCA2"), colnames(exprFinal_mLDA_cross))
try2_new  <- mLDA.pair(Xnew, caseControlResponse, Xnew, Z=NULL, Z.new=NULL, pair=c(1,0), strong.X.set=c(1:50,927,625), tau=52, alpha=0.8, nu=200, d=2, nb=5)   
saveRDS(try2_new,file=paste("case_contral_mLDA_V2", ".rds", sep=""))
pred <- as.numeric(try2_new$PredClass)
length(which(pred!=Y))
try2_new <- readRDS(file=paste("case_contral_mLDA_V2", ".rds", sep=""))
#write.csv(as.data.frame(try2_new[["network.nodes"]]),"case_contral_mLDA_V2.csv")
#write.csv(as.data.frame(try2_new[["MIset"]]),"case_contral_mLDA_V2_strong.csv")
survival_mLDA_network <- NULL
for (i in 1:22) {
  temp <- as.data.frame(try2_new[["local.network_list"]][[i]])
  temp$gene <- rownames(temp)
  temp$net <- paste("Network",i, sep = "")
  temp <- temp[,-1]
  survival_mLDA_network <- rbind(survival_mLDA_network,temp)
}
write.csv(survival_mLDA_network,"survival_mLDA_network.csv")

################################################################
################################################################
#case-contral study using the turnning parameters decided by CV
#cross validation
caseControlResponse <- rep(0, (dim(exprFinal)[2]-2))
caseControlResponse[which(sample_vec == "TCGA")] <- 1
try2  <- mLDA.cv(exprFinal_mLDA_cross, caseControlResponse, Z=NULL, fold=5, seed=1, tau=50, alpha=0.9, nu=150, d=2, nb=5)
saveRDS(try2,file=paste("case_contral_crossvalidation_mLDA", ".rds", sep=""))
#try2 = readRDS(file=paste("case_contral_crossvalidation_mLDA", ".rds", sep=""))

glm.pair.crossvalidation <-
  function(X, Y, sepn){
    valid_output = NULL
    
    for (i in 1:sepn) {
      n = dim(X)[2]
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
      
      #### glm/glmnet fit ####
      library(tidyverse)
      library(glmnet)
      #coef_fit <- glm(y_Train~X_train,family = binomial) %>% coef
      #temp <- glmnet(x=X_train,y=y_Train, alpha = 1,lambda=0.01, family = binomial(link = "logit"), relax=TRUE) #lasso
      temp <- glmnet(x=X_train,y=y_Train, alpha = 0,lambda=0.01, family = binomial(link = "logit"), relax=TRUE) #ridge
      #temp <- glmnet(x=X_train,y=y_Train, alpha = 0.5,lambda=0.01, family = binomial(link = "logit"), relax=TRUE) #elastic net
      coef_fit <- rbind(temp[["a0"]],as.matrix(temp[["beta"]]))
      #length(which(temp[["beta"]]!=0))
      #p_prediction <- exp(cbind(1,X_test)%*%coef_fit)/(1+exp(cbind(1,X_test)%*%coef_fit))
      p_prediction <- exp(cbind(1,X_test)%*%coef_fit)/(1+exp(cbind(1,X_test)%*%coef_fit))
      y_prediction <- rbinom(dim(X_test)[1], 1, p_prediction)
      temp <- length(which(y_prediction != y_Test))
      valid_output <- rbind(valid_output,temp)
    }
    return(valid_output)
  }
result <-glm.pair.crossvalidation(X=t(Xnew[,1:500]), Y=caseControlResponse, sepn=5)
result
result <-glm.pair.crossvalidation(X=t(Xnew), Y=caseControlResponse, sepn=5)
result

################################################################
################################################################
#Permutation pvalue Caculation
setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control")
try2 = readRDS(file=paste("case_contral_mLDA_V2", ".rds", sep=""))
All_CpG = rownames(as.data.frame(try2[["network.nodes"]]))
X_All <- exprFinal_mLDA_cross[,colnames(exprFinal_mLDA_cross)%in%All_CpG]
X_All <- exprFinal_mLDA_cross
#mLDA.permutation.test(X, Y, Z=NULL, K=2, select.set, Iff, Omega.X, permN=100000)
# scale dataset 
# for this dataset exprFinal_mLDA_cross has been normalized
#X_All = t(scale(t(X_All), center = TRUE, scale = TRUE))
All_y = as.data.frame(caseControlResponse)
X_All = cbind(All_y, X_All)
#X_All=as.data.frame(t(X_All))
#Omega.r <- try2[["network.Omega"]]
Omega.r <- my.inv(cov(X_All[,-1]))

ISresample_12M <- NULL
#Permutation_times = 100000
for (i in 1:100000){
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

ISresample_12M_A <- matrix(rep(abs(ISresample_12M[,1]), 100000), dim(ISresample_12M)[1], 100000) 
#Here abs function is making sure there is no error for further steps
ISresample_12M_B <- apply(ISresample_12M[, -1], 1, abs)
ISresample_12M_B <- apply(ISresample_12M[, -1], 2, abs)
CCC <- ISresample_12M_B-ISresample_12M_A 
ddd <- apply(CCC, 1, function(x) length(which(x>0)))/100000 
pvalue_CpG_sites[,2] = ddd
pvalue_CpG_sites = as.data.frame(pvalue_CpG_sites)
pvalue_CpG_sites$V2[pvalue_CpG_sites$V2==0] <- 10e-10
write.csv(pvalue_CpG_sites, "C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/case_control_mLDA_sig.csv")

######################################################################
###########Hotelling’s test on a selected network###################
######################################################################
nodes_name <- rownames(as.data.frame(try2_new[["local.network_list"]][[22]]))
temp <- t(ExprDD)[,rownames(ExprDD)%in%nodes_name]
#rownames(temp)<-1:dim(temp)[1]   #only in hotelling
XpairSel <- temp

ypair <- caseControlResponse
ypair[which(caseControlResponse==0)] <- 1
ypair[which(caseControlResponse==1)] <- 2

temp1 <- cbind(ypair,XpairSel)
t.test(XpairSel~ypair, data=temp1)
temp1_g1 <- as.data.frame(subset(temp1, ypair==1))
temp1_g1 <- temp1_g1[,-1]
temp1_g1_expr <- as.data.frame(apply(temp1_g1,2,sum))   #only in hotelling

temp1_g2 <- as.data.frame(subset(temp1, ypair==2))
temp1_g2 <- temp1_g2[,-1]
temp1_g2_expr <- as.data.frame(apply(temp1_g2,2,mean))   #only in hotelling

temp1 <- cbind(temp1_g1_expr,temp1_g2_expr)
colnames(temp1)<-c("GTEx","TCGA")
muH0 <- c(0, 0)
library(DescTools)
HotellingsT2Test(x=temp1, mu=muH0)   ##only in hotelling



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

case_control_edgeR_GSEA_KEGG_output <- read_excel("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/GSEA_CC/case_control_edgeR_GSEA_KEGG_output.xlsx", 
                                                  +     sheet = "Sheet1")
colnames(case_control_edgeR_GSEA_KEGG_output)
case_control_edgeR_GSEA_KEGG_output <- case_control_edgeR_GSEA_KEGG_output[,c(2,4,5,7)]
colnames(case_control_edgeR_GSEA_KEGG_output) <-c("GS","SIZE","ES","NOM_pval")
ggplot(case_control_edgeR_GSEA_KEGG_output) +
  geom_point(aes(x = ES, y = reorder(GS, ES), color = NOM_pval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")+
  theme(text = element_text(size = 12))

#mLDA output
case_control_mLDA_GSEA_output_upregulated <- read_excel("GSEA_CC/case_control_mLDA_GSEA_output_upregulated.xlsx")
colnames(case_control_edgeR_GSEA_KEGG_output )
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

###################V2 of GSEA
##### GO ontology
case_control_edgeR_GSEA_GO_output <- read_excel("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/GSEA_CC/case_control_edgeR_GSEA_GO_output.xlsx")
#GSEAdata = subset(case_control_edgeR_GSEA_GO_output, SIZE!=1)
#GSEAdata = subset(case_control_edgeR_GSEA_GO_output, ES>0)
#GSEAdata = subset(case_control_edgeR_GSEA_GO_output, ES<0)
colnames(GSEAdata)[7] = "NOM_pval"
colnames(GSEAdata)[8] = "FDR_qval"
ggplot(GSEAdata) +
  geom_point(aes(x = ES, y = reorder(NAME, ES), color = FDR_qval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")

##### KEGG
case_control_edgeR_GSEA_KEGG_output <- read_excel("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Results/Case_control/GSEA_CC/case_control_edgeR_GSEA_KEGG_output.xlsx")
#GSEAdata = subset(case_control_edgeR_GSEA_KEGG_output, SIZE!=1)
#GSEAdata = subset(case_control_edgeR_GSEA_KEGG_output, ES>0)
#GSEAdata = subset(case_control_edgeR_GSEA_KEGG_output, ES<0)
colnames(GSEAdata)[7] = "NOM_pval"
colnames(GSEAdata)[8] = "FDR_qval"
ggplot(GSEAdata) +
  geom_point(aes(x = ES, y = reorder(NAME, ES), color = FDR_qval, size=SIZE)) +
  scale_color_gradient(low ="green",high ="red")



######################################################################
###################   Visualization of gene network ##
######################################################################
nodes_name <- rownames(as.data.frame(try2_new[["local.network_list"]][[18]]))
temp <- exprFinal_mLDA_cross[,colnames(exprFinal_mLDA_cross)%in%nodes_name]
rownames(temp)<-1:dim(temp)[1]
XpairSel <- temp
library(igraph)
g <- graph.adjacency(
  as.matrix(as.dist(cor(temp, method="pearson"))),
  mode="undirected",
  weighted=TRUE,
  diag=FALSE
)

# Simplfy the adjacency object
g <- simplify(g, remove.multiple=TRUE, remove.loops=TRUE)

# Colour negative correlation edges as blue
E(g)[which(E(g)$weight<0)]$color <- "darkblue"

# Colour positive correlation edges as red
E(g)[which(E(g)$weight>0)]$color <- "darkred"

# Convert edge weights to absolute values
E(g)$weight <- abs(E(g)$weight)

# Remove edges below absolute Pearson correlation 0.8
g <- delete_edges(g, E(g)[which(E(g)$weight<0.5)])

# Remove any vertices remaining that have no edges
g <- delete.vertices(g, degree(g)==0)

# Assign names to the graph vertices (optional)
V(g)$name <- V(g)$name

# Change shape of graph vertices
V(g)$shape <- "sphere"

# Change colour of graph vertices
V(g)$color <- "skyblue"

# Change colour of vertex frames
V(g)$vertex.frame.color <- "white"

# Scale the size of the vertices to be proportional to the level of expression of each gene represented by each vertex
# Multiply scaled vales by a factor of 10
vSizes <- (temp1 + 1.0) * 10

# Amplify or decrease the width of the edges
edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
mst <- mst(g, algorithm="prim")

# Plot the tree object
plot(
  mst,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  vertex.size=vSizes,
  vertex.label.dist=-0.5,
  vertex.label.color="black",
  asp=FALSE,
  vertex.label.cex=0.6,
  edge.width=edgeweights,
  edge.arrow.mode=0,
  main="Network #18")
