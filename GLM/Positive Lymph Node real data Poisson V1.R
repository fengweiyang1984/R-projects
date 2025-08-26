##################################################################################################
##### case-control classification analysis
##################################################################################################
##################################################################################################
setwd("C:/Users/z021w783/Desktop/GRA/Disertation/Real data example")
#### single gene based DE analysis
#### 
sample_vec <- sapply(colnames(exprFinal), function(x) unlist(strsplit(x, "-"))[[1]][1]) 
table(sample_vec)
#sample_vec
#gene GTEX   id TCGA 
#1  178    1 1092 
rownames(exprFinal) <- exprFinal[,2]
BRCA_clinicalMatrix <- read.delim("TCGA.BRCA.sampleMap/BRCA_clinicalMatrix")

##########Extract the known genes in disease related KEGG pathways
library("KEGGREST")
namestrue <- as.data.frame(NULL)
#Get the list of numbers, gene symbols and gene description
names <- keggGet("hsa05224")[[1]]$GENE    #BC pathway
namesodd <-  names[seq(0,length(names),2)]
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

#names <- keggGet("hsa05200")[[1]]$GENE   #pathway in cancer
#namesodd <-  names[seq(0,length(names),2)]
#temp <- as.data.frame(gsub("\\;.*","",namesodd))
#namestrue <- rbind(namestrue,temp)

temp <- toupper(namestrue[,1])
#subset X matrix based on gene symbols from namestrue
X_input <- t(exprFinal[,-1:-2])
X_input <- apply(X_input,2,scale)
X_input1 <- X_input[,colnames(X_input)%in%temp]
dim(X_input1)
rownames(X_input1) <- colnames(exprFinal[,-1:-2])
X_input1 <- X_input1[rownames(X_input1)%in%BRCA_clinicalMatrix$sampleID,]
dim(X_input1)
y_input <- BRCA_clinicalMatrix$number_of_lymphnodes_positive_by_he[which(BRCA_clinicalMatrix$sampleID%in%rownames(X_input1))]
y_input[which(is.na(y_input))] <- 0
length(y_input)

#####cross validation to determine the best lambda in group lasso
IRLS.IP.gglasso.crossvalidation <-
  function(X2, y2, n, sepn =5){
    valid_result <- NULL
    for (k in seq(0.1, 3, length.out = 30)) {
    valid_output = NULL
    for (tt in 1:sepn) {
      #### Preparing the training & testing dataset
      lowb1 = ceiling(n/sepn*(tt-1))+1
      upb1 = min(ceiling(n/sepn*tt),n)
      data_set1 = X2
      
      X_train = data_set1[-lowb1:-upb1,]
      X_test = data_set1[lowb1:upb1,]
      
      y_Train = y2[-lowb1:-upb1]
      y_Test =  y2[lowb1:upb1]
      
      XtY <- cor(X_train,y_Train)
      XtY[which(abs(XtY)<0.4)]<-0
      XtY <- as.data.frame(XtY)
      XtX <- cor(X_train)
      Sigma <- cov(X_train)
      Sigma <- round(Sigma,digits = 5)
      library(matlib)
      Sigma_inv <- solve(Sigma)
      dim(Sigma_inv)
      Sigma_inv[which(XtX<=0.5)] <- 0

      glm_mine <- IRLS_IP_gglasso(X=X_train,y=y_Train,lambda_lasso=k)
      temp.lasso <- glmnet(x=X_train,y=y_Train, alpha = 0.5,lambda=k, poisson(link = "log"),intercept = T)
      p_prediction <- predict(temp.lasso,  newx=X_test, type="response")
      
      p_prediction <- X_test%*%glm_mine
      p_prediction[which(p_prediction<=0)] <- 0
      y_prediction <- rpois(dim(X_test)[1],p_prediction)
      temp_result1 <- mean((y_Test-y_prediction)^2)/(dim(X_test)[1]-length(which(glm_mine!=0))) #MSE  mean estimation sum of squares
      temp_result2 <- 1-(sum((y_prediction-y_Test)^2)/(dim(X_test)[1]-length(which(glm_mine!=0))))/(sum((y_Test-mean(y_Test))^2)/(dim(X_test)[1]-1)) #adjusted R2
      
      valid_output <- rbind(valid_output,cbind(temp_result1,temp_result2))
    }
  temp1 <- as.vector(c(k,apply(valid_output,2,mean)))
  valid_result <- rbind(valid_result,temp1)
  print(k)
  }
  return(valid_result)
}
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/GLM Iterative Proximal (IP) Algorithm Sparse Regression for Poisson V_gglasso5.R")
p=dim(X_input1)[2]
library("glmnet")
valid_result <- IRLS.IP.gglasso.crossvalidation(X2=X_input1, y2=y_input, n=nrow(X_input1),sepn =5)
write.csv(valid_result,"valid_result.csv")

#correlation betweeen X and y
XtY <- cor(X_input1,y_input)
XtY[which(abs(XtY)<0.4)]<-0
XtY <- as.data.frame(XtY)

#correlation of X input
XtX <- cor(X_input1)
#covariance of X
Sigma <- cov(X_input1)
Sigma <- round(Sigma,digits = 5)
library(matlib)
#precision matrix
Sigma_inv <- solve(Sigma)
dim(Sigma_inv)
#use a threshold on correlation of X to threshold precision matrix
Sigma_inv[which(XtX<=0.5)] <- 0

glm_mine <- IRLS_IP_gglasso(X=X_input1,y=y_input,lambda_lasso=0.3)
length(which(glm_mine!=0))
rownames(glm_mine) <- colnames(X_input1)
write.csv(glm_mine,"valid_result.csv")

temp.lasso <- glmnet(x=X_input1,y=y_input, alpha = 0.5,lambda=0.4, poisson(link = "log"),intercept = F)
glm_laaso <- as.matrix(temp.lasso[["beta"]])
length(which(glm_laaso!=0))
write.csv(glm_laaso,"valid_result.csv")

###################   Visualization of pathways
library(pathview)

#rownames are hsa id for genes, first column is prerank score.
pathway_comparison <- as.matrix(Poisson_result)

pv.out <- pathview(gene.data = pathway_comparison[,1], pathway.id = "05224",species = "hsa", 
                   out.suffix = "gse16873.2layer", kegg.native = T,same.layer = F,sign.pos = "bottomleft")

