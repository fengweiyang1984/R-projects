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
sample_vec<-sample_vec[-1:-2]
caseControlResponse <- rep(0, (dim(exprFinal)[2]-2))
caseControlResponse[which(sample_vec == "TCGA")] <- 1
rownames(exprFinal) <- exprFinal[,2]

##########Extract the known genes in disease related KEGG pathways
library("KEGGREST")
namestrue <- as.data.frame(NULL)
#Get the list of numbers, gene symbols and gene description
names <- keggGet("mmu05224")[[1]]$GENE
#Delete the gene number by deleting every other line
namesodd <-  names[seq(0,length(names),2)]
#Create a substring deleting everything after the ; on each line (this deletes the gene description).
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

names <- keggGet("mmu05206")[[1]]$GENE
namesodd <-  names[seq(0,length(names),2)]
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

names <- keggGet("mmu04151")[[1]]$GENE
namesodd <-  names[seq(0,length(names),2)]
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

names <- keggGet("hsa05212")[[1]]$GENE
namesodd <-  names[seq(0,length(names),2)]
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

names <- keggGet("hsa05200")[[1]]$GENE
namesodd <-  names[seq(0,length(names),2)]
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

temp <- toupper(namestrue[,1])
#subset X matrix based on gene symbols from namestrue
X_input <- t(exprFinal[,-1:-2])
X_input <- apply(X_input,2,scale)
X_input <- X_input[,colnames(X_input)%in%temp]
dim(X_input)
y_input <- as.integer(caseControlResponse)

#####cross validation to determine the best lambda in group lasso
IRLS.IP.gglasso.crossvalidation <-
  function(X2, y2, n, sepn =5,lambda = 5){
    valid_result <- NULL
    for (k in seq(0.21, 0.3, length.out = 9)) {
    valid_output = NULL
    for (tt in 1:sepn) {
      #### Preparing the training & testing dataset
      lowb1 = ceiling(178/sepn*(tt-1))+1
      upb1 = min(ceiling(178/sepn*tt),178)
      lowb2 = ceiling((n-178)/sepn*(tt-1))+1
      upb2 = min(ceiling((n-178)/sepn*tt),n-180)
      data_set1 = X2[1:178,]
      data_set2 = X2[-1:-178,]
      
      data_set1_Train = data_set1[-lowb1:-upb1,]
      data_set2_Train = data_set2[-lowb2:-upb2,]
      X_train = rbind(data_set1_Train,data_set2_Train)
      data_set1_Test = data_set1[lowb1:upb1,]
      data_set2_Test = data_set2[lowb2:upb2,]
      X_test = rbind(data_set1_Test,data_set2_Test) 
      response_set1 = y2[1:178]
      response_set2 = y2[-1:-178]
      
      y_Train = as.vector(c(response_set1[-lowb1:-upb1],response_set2[-lowb2:-upb2]))
      y_Test =  as.vector(c(response_set1[lowb1:upb1],response_set2[lowb2:upb2]))
      
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
      X<-X_train
      y<-y_Train
      glm_mine <- IRLS_IP_gglasso(X=X,y=y,lambda_lasso=k)
      p_prediction <- exp(X_test%*%glm_mine)/(1+exp(X_test%*%glm_mine))

      y_prediction <- rep(2,dim(X_test)[1])
      y_prediction[which(p_prediction>=0.5)] <- 1
      y_prediction[which(p_prediction<=0.5)] <- 0
      
      temp <- length(which(y_prediction==y_Test))/length(y_Test)
      valid_output <- rbind(valid_output,temp)
    }
  temp1 <- as.vector(c(k,mean(valid_output)))
  valid_result <- rbind(valid_result,temp1)
  print(k)
  }
  return(valid_result)
}
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/GLM Iterative Proximal (IP) Algorithm Sparse Regression for logistic V_gglasso5.R")
valid_result <- IRLS.IP.gglasso.crossvalidation(X2=X_input, y2=y_input, n=ncol(X_input),sepn =5)
write.csv(valid_result,"valid_result.csv")

Crossvalidation_results$lambda1 = 1-Crossvalidation_results$lambda
ggplot(data = Crossvalidation_results, aes(x = lambda1, y = Prediction_accuracy)) +
  geom_line()

#correlation betweeen X and y
XtY <- cor(X_input,y_input)
XtY[which(abs(XtY)<0.4)]<-0
XtY <- as.data.frame(XtY)

#correlation of X input
XtX <- cor(X_input)
#covariance of X
Sigma <- cov(X_input)
Sigma <- round(Sigma,digits = 5)
library(matlib)
#precision matrix
Sigma_inv <- solve(Sigma)
dim(Sigma_inv)
#use a threshold on correlation of X to threshold precision matrix
Sigma_inv[which(XtX<=0.5)] <- 0

#model fitting
#group lasso
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/GLM Iterative Proximal (IP) Algorithm Sparse Regression for logistic V_gglasso5.R")
try(glm_mine <- IRLS_IP_gglasso(X=X_input,y=y_input,lambda_lasso=0.03))

#sparse group lasso
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/GLM Iterative Proximal (IP) Algorithm Sparse Regression for logistic V_sparsegl5.R")
try(glm_mine1 <- IRLS_IP_sparsegl(X=X_input,y=y_input,lambda_lasso=0))

library("glmnet") 
#lasso
temp <- glmnet(x=X_input,y=y_input, alpha = 1,lambda=0.05, poisson(link = "log"), relax=T,intercept = F)
glm_laaso <- as.matrix(temp[["beta"]])
glm_laaso_bo <- temp[["a0"]]

#elastic net
temp <- glmnet(x=X_input,y=y_input, alpha = 0.5,lambda=0.1, poisson(link = "log"), relax=T,intercept = F)
glm_EN <- as.matrix(temp[["beta"]])
glm_EN_bo <- temp[["a0"]]

###
temp <- cbind(glm_mine,glm_mine1,glm_laaso,glm_EN) 
colnames(temp) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_EN")
beta_estimation <- temp
y_prediction <- NULL
for (i in 1:ncol(beta_estimation)){
  temp <- exp(X_input%*%beta_estimation[,i])/(1+exp(X_input%*%beta_estimation[,i]))
  y_prediction <- cbind(y_prediction,temp)
}
colnames(y_prediction) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_EN")
y_prediction <- as.data.frame(y_prediction)
y_prediction$true_outcome <- caseControlResponse


# true positive (hit) rate
tpr <- function(pred, actual) {
  res <- data.frame(pred, actual)
  sum(res$actual == 1 & res$pred == 1) / sum(actual == 1)
}

# false positive rate
fpr <- function(pred, actual) {
  res <- data.frame(pred, actual)
  sum(res$actual == 0 & res$pred == 1) / sum(actual == 0)
}

#################use different thresholding to calculate tpr and fpr
#################they would be adopted in ROC curve for logistic
###IRLS_IP_gglasso
temp <- matrix(c(0, 0, 0, 1, 1, 1), ncol=3, byrow=TRUE)
colnames(temp) <- c('thresholding','tpr','fpr')
temp <- as.data.frame(temp)
for (i in 1:90){
  thresholding = i/1800
  thresholding1 = i/300
  temp1 <- y_prediction_summary
  temp1$IRLS_IP_gglasso_pred <- rep(2,nrow(temp1))
  temp1$IRLS_IP_gglasso_pred[which(temp1$IRLS_IP_gglasso>=(0.8-thresholding1))]<- 1
  tpr = length(which(temp1$IRLS_IP_gglasso_pred[which(temp1$true_outcome==1)]==1)) / length(which(temp1$true_outcome == 1))
  temp1$IRLS_IP_gglasso_pred[which(temp1$IRLS_IP_gglasso<=(0.6))]<- 1
  temp1$IRLS_IP_gglasso_pred[which(temp1$IRLS_IP_gglasso<=(0.05-thresholding))]<- 0
  fpr = length(which(temp1$IRLS_IP_gglasso_pred[which(temp1$true_outcome==0)]==1)) / length(which(temp1$true_outcome == 0))
  temp <- rbind(temp,c(thresholding,tpr,fpr))
}
write.csv(temp,"temp.csv")

temp <- matrix(c(0, 0, 0, 1, 1, 1), ncol=3, byrow=TRUE)
colnames(temp) <- c('thresholding','tpr','fpr')
temp <- as.data.frame(temp)
for (i in 1:90){
  thresholding = i/1800
  thresholding1 = i/300
  temp1 <- y_prediction_summary
  temp1$IRLS_IP_sparsegl_pred <- rep(2,nrow(temp1))
  temp1$IRLS_IP_sparsegl_pred[which(temp1$IRLS_IP_sparsegl>=(0.8-thresholding1))]<- 1
  tpr = length(which(temp1$IRLS_IP_sparsegl_pred[which(temp1$true_outcome==1)]==1)) / length(which(temp1$true_outcome == 1))
  temp1$IRLS_IP_sparsegl_pred[which(temp1$IRLS_IP_sparsegl<=(0.6))]<- 1
  temp1$IRLS_IP_sparsegl_pred[which(temp1$IRLS_IP_sparsegl<=(0.05-thresholding))]<- 0
  fpr = length(which(temp1$IRLS_IP_sparsegl_pred[which(temp1$true_outcome==0)]==1)) / length(which(temp1$true_outcome == 0))
  temp <- rbind(temp,c(thresholding,tpr,fpr))
}
write.csv(temp,"temp.csv")

















library(plotROC)
ggplot(temp, aes(x = fpr, y = tpr, color = method)) + geom_line()
ggplot(temp, aes(x = fpr, y = tpr)) + geom_line()
