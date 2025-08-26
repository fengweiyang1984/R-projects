#source("https://neuroconductor.org/neurocLite.R")
#neuro_install("neurobase", release = "stable")
setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 3/Results")
Logistic_Data <- regionM[1:184,]   #AD_image 53    MCI_image  131


cor_Logistic_Data <- cor(Logistic_Data)
#cor_Logistic_Data[which(abs(cor_Logistic_Data)<0.8, arr.ind=TRUE)] <- 0
cor_Logistic_Data0 <- cor_Logistic_Data
cor_Logistic_Data0[which(pre_network==0)]<-0

X_input <- Logistic_Data
X_input <- apply(X_input,2,scale)
dim(X_input)
caseControlResponse <- rep(0, dim(X_input)[1])
caseControlResponse[1:53] <- 1
y_input <- as.integer(caseControlResponse)

################################################################
################################################################
################cross validation for lambda#####################
################################################################
################################################################

##gglasso
IRLS.IP.gglasso.crossvalidation <-
  function(X2, y2, n, sepn =10,lambda = 5){
    valid_result <- NULL
    for (k in seq(0.01, 0.15, length.out = 15)) {
      valid_output = NULL
      for (tt in 2:sepn) {
        #### Preparing the training & testing dataset
        lowb1 = ceiling(53/sepn*(tt-1))+1
        upb1 = min(ceiling(53/sepn*tt),53)
        lowb2 = ceiling((n-53)/sepn*(tt-1))+1
        upb2 = min(ceiling((n-53)/sepn*tt),n-180)
        data_set1 = X2[1:53,]
        data_set2 = X2[-1:-53,]
        
        data_set1_Train = data_set1[-lowb1:-upb1,]
        data_set2_Train = data_set2[-lowb2:-upb2,]
        X_train = rbind(data_set1_Train,data_set2_Train)
        data_set1_Test = data_set1[lowb1:upb1,]
        data_set2_Test = data_set2[lowb2:upb2,]
        X_test = rbind(data_set1_Test,data_set2_Test) 
        response_set1 = y2[1:53]
        response_set2 = y2[-1:-53]
        
        y_Train = as.vector(c(response_set1[-lowb1:-upb1],response_set2[-lowb2:-upb2]))
        y_Test =  as.vector(c(response_set1[lowb1:upb1],response_set2[lowb2:upb2]))
        
        XtY <- cor(X_train,y_Train)
        XtY[which(abs(XtY)<0.4)]<-0
        XtY <- as.data.frame(XtY)
        XtX <- pre_network1
        Sigma <- cov(X_train)
        Sigma <- round(Sigma,digits = 5)
        library(matlib)
        Sigma_inv <- solve(Sigma)
        dim(Sigma_inv)
        Sigma_inv[which(abs(XtX)==0, arr.ind=TRUE)] <- 0
        #Sigma_inv[which(pre_network==0)]<-0
        X<-X_train
        y<-y_Train
        
        temp <- glmnet(x=X,y=y, alpha = 1,lambda=k, poisson(link = "log"))
        glm_laaso <- as.matrix(temp[["beta"]])
        intercept <- temp[["a0"]]
        
        glm_mine <- IRLS_IP_gglasso(X=X,y=y,lambda_lasso=k)
        #p_prediction <- exp(X_test%*%glm_mine)/(1+exp(X_test%*%glm_mine))
        p_prediction <- exp(X_test%*%glm_mine+as.matrix(rep(intercept,nrow(X_test))))/(1+exp(X_test%*%glm_mine+as.matrix(rep(intercept,nrow(X_test)))))
        
        y_prediction <- rep(0,dim(X_test)[1])
        y_prediction[which(p_prediction>=0.5)] <- 1
        
        temp <- length(which(y_prediction==y_Test))/length(y_Test)
        valid_output <- rbind(valid_output,temp)
      }
      temp1 <- as.vector(c(k,mean(valid_output)))
      valid_result <- rbind(valid_result,temp1)
      print(k)
    }
    return(valid_result)
  }
library(glmnet)
source("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 3/Code/GLM Iterative Proximal (IP) Algorithm Sparse Regression for logistic V_gglasso5.R")
valid_result <- IRLS.IP.gglasso.crossvalidation(X2=X_input, y2=y_input, n=nrow(X_input),sepn =5)
write.csv(valid_result,"valid_result.csv")


#lasso
valid_result <- NULL
for (k in seq(0.01, 0.15, length.out = 15)) {
  valid_output = NULL
  for (tt in 1:sepn) {
    #### Preparing the training & testing dataset
    lowb1 = ceiling(53/sepn*(tt-1))+1
    upb1 = min(ceiling(53/sepn*tt),53)
    lowb2 = ceiling((n-53)/sepn*(tt-1))+1
    upb2 = min(ceiling((n-53)/sepn*tt),n-180)
    data_set1 = X2[1:53,]
    data_set2 = X2[-1:-53,]
    
    data_set1_Train = data_set1[-lowb1:-upb1,]
    data_set2_Train = data_set2[-lowb2:-upb2,]
    X_train = rbind(data_set1_Train,data_set2_Train)
    data_set1_Test = data_set1[lowb1:upb1,]
    data_set2_Test = data_set2[lowb2:upb2,]
    X_test = rbind(data_set1_Test,data_set2_Test) 
    response_set1 = y2[1:53]
    response_set2 = y2[-1:-53]
    
    y_Train = as.vector(c(response_set1[-lowb1:-upb1],response_set2[-lowb2:-upb2]))
    y_Test =  as.vector(c(response_set1[lowb1:upb1],response_set2[lowb2:upb2]))
    X<-X_train
    y<-y_Train
    
    fit <- glmnet(x=X_input,y=y_input, alpha = 1,lambda=k, family ="binomial")
    p_prediction <- exp(X_test%*%fit$beta)/(1+exp(X_test%*%fit$beta))
    
    y_prediction <- rep(0,dim(X_test)[1])
    y_prediction[which(p_prediction>=0.5)] <- 1
    
    temp <- length(which(y_prediction==y_Test))/length(y_Test)
    valid_output <- rbind(valid_output,temp)
  }
  temp1 <- as.vector(c(k,mean(valid_output)))
  valid_result <- rbind(valid_result,temp1)
  print(k)
}
write.csv(valid_result,"valid_result.csv")

#elastic net
valid_result <- NULL
for (k in seq(0.01, 0.15, length.out = 15)) {
  valid_output = NULL
  for (tt in 1:sepn) {
    #### Preparing the training & testing dataset
    lowb1 = ceiling(53/sepn*(tt-1))+1
    upb1 = min(ceiling(53/sepn*tt),53)
    lowb2 = ceiling((n-53)/sepn*(tt-1))+1
    upb2 = min(ceiling((n-53)/sepn*tt),n-180)
    data_set1 = X2[1:53,]
    data_set2 = X2[-1:-53,]
    
    data_set1_Train = data_set1[-lowb1:-upb1,]
    data_set2_Train = data_set2[-lowb2:-upb2,]
    X_train = rbind(data_set1_Train,data_set2_Train)
    data_set1_Test = data_set1[lowb1:upb1,]
    data_set2_Test = data_set2[lowb2:upb2,]
    X_test = rbind(data_set1_Test,data_set2_Test) 
    response_set1 = y2[1:53]
    response_set2 = y2[-1:-53]
    
    y_Train = as.vector(c(response_set1[-lowb1:-upb1],response_set2[-lowb2:-upb2]))
    y_Test =  as.vector(c(response_set1[lowb1:upb1],response_set2[lowb2:upb2]))
    X<-X_train
    y<-y_Train
    
    fit <- glmnet(x=X_input,y=y_input, alpha = 0.5,lambda=k, family ="binomial")
    p_prediction <- exp(X_test%*%fit$beta)/(1+exp(X_test%*%fit$beta))
    
    y_prediction <- rep(0,dim(X_test)[1])
    y_prediction[which(p_prediction>=0.5)] <- 1
    
    temp <- length(which(y_prediction==y_Test))/length(y_Test)
    valid_output <- rbind(valid_output,temp)
  }
  temp1 <- as.vector(c(k,mean(valid_output)))
  valid_result <- rbind(valid_result,temp1)
  print(k)
}
write.csv(valid_result,"valid_result.csv")

#ridge
valid_result <- NULL
for (k in seq(0.01, 0.15, length.out = 15)) {
  valid_output = NULL
  for (tt in 1:sepn) {
    #### Preparing the training & testing dataset
    lowb1 = ceiling(53/sepn*(tt-1))+1
    upb1 = min(ceiling(53/sepn*tt),53)
    lowb2 = ceiling((n-53)/sepn*(tt-1))+1
    upb2 = min(ceiling((n-53)/sepn*tt),n-180)
    data_set1 = X2[1:53,]
    data_set2 = X2[-1:-53,]
    
    data_set1_Train = data_set1[-lowb1:-upb1,]
    data_set2_Train = data_set2[-lowb2:-upb2,]
    X_train = rbind(data_set1_Train,data_set2_Train)
    data_set1_Test = data_set1[lowb1:upb1,]
    data_set2_Test = data_set2[lowb2:upb2,]
    X_test = rbind(data_set1_Test,data_set2_Test) 
    response_set1 = y2[1:53]
    response_set2 = y2[-1:-53]
    
    y_Train = as.vector(c(response_set1[-lowb1:-upb1],response_set2[-lowb2:-upb2]))
    y_Test =  as.vector(c(response_set1[lowb1:upb1],response_set2[lowb2:upb2]))
    X<-X_train
    y<-y_Train
    
    fit <- glmnet(x=X_input,y=y_input, alpha = 0,lambda=k, family ="binomial")
    p_prediction <- exp(X_test%*%fit$beta)/(1+exp(X_test%*%fit$beta))
    
    y_prediction <- rep(0,dim(X_test)[1])
    y_prediction[which(p_prediction>=0.5)] <- 1
    
    temp <- length(which(y_prediction==y_Test))/length(y_Test)
    valid_output <- rbind(valid_output,temp)
  }
  temp1 <- as.vector(c(k,mean(valid_output)))
  valid_result <- rbind(valid_result,temp1)
  print(k)
}
write.csv(valid_result,"valid_result.csv")

set.seed(10101)
cvob1 = cv.glmnet(x=X_input, y=y_input, family = "binomial", type.measure = "class", nfolds = 5, alpha = 1)
plot(cvob1)
cvob1$lambda.1se
cvob1$lambda.min
write.csv(cbind(cvob1$lambda,cvob1$cvm), file = "cvlasso.csv")
cvob2 = cv.glmnet(x=X_input, y=y_input, family = "binomial", type.measure = "class", nfolds = 5, alpha = 0)
plot(cvob2)
cvob2$lambda.min
write.csv(cbind(cvob2$lambda,cvob2$cvm), file = "cvridge.csv")
cvob3 = cv.glmnet(x=X_input, y=y_input, family = "binomial", type.measure = "class", nfolds = 5, alpha = 0.5)
plot(cvob3)
cvob3$lambda.min
write.csv(cbind(cvob3$lambda,cvob3$cvm), file = "cvEN.csv")

################################################################
################################################################
################cross validation for lambda#####################
################################################################
################################################################

XtY <- cor(X_input,y_input)
XtY[which(abs(XtY)<0.6)]<-0
XtY <- as.data.frame(XtY)
XtX <- cor(X_input)
Sigma <- cov(X_input)
Sigma <- round(Sigma,digits = 5)
library(matlib)
Sigma_inv <- solve(Sigma)
dim(Sigma_inv)
#Sigma_inv[which(abs(XtX)<0.6, arr.ind=TRUE)] <- 0
Sigma_inv[which(pre_network1==0)]<-0
glm_mine <- IRLS_IP_gglasso(X=X_input,y=y_input,lambda_lasso=0.07)
length(which(glm_mine!=0))
region_code[which(glm_mine!=0),]
glm_mine <- as.data.frame(glm_mine)
rownames(glm_mine) <- region_code[,3]
write.csv(glm_mine,"HDnetGLM.csv")

glm_lasso = glmnet(x=X_input, y=y_input, family = "binomial", alpha = 1, lambda=0.062)
glm_laaso <- as.matrix(glm_lasso[["beta"]])
length(which(glm_laaso!=0))
write.csv(region_code[which(glm_laaso!=0),],"Lasso.csv")
glm_ridge = glmnet(x=X_input, y=y_input, family = "binomial", alpha = 0, lambda=0.175)
glm_ridge <- as.matrix(glm_ridge[["beta"]])
length(which(glm_ridge!=0))
glm_EN = glmnet(x=X_input, y=y_input, family = "binomial", alpha = 0.5, lambda=0.093)
glm_EN <- as.matrix(glm_EN[["beta"]])
length(which(glm_EN!=0))

glm_mine_nonzero <- subset(glm_mine,V1!=0)
X_input1 <- X_input[,colnames(X_input)%in%region_code[which(glm_mine!=0),2]]
X_input1 <- as.data.frame(cbind(X_input1,y_input))
colnames(X_input1)
pvalue_glm_mine <- glm(y_input~`4001`,data =X_input1, family = "binomial")
summary(pvalue_glm_mine)
