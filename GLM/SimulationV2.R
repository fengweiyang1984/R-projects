##########################################################################
############################# Simulation #################################
##########################################################################

library(tidyverse)
library(MASS)
#p is # of variables
n=500
p=500
sep <- 50

Sigma_sub_inv <- matrix(0,p/sep,p/sep)   
for(i in 1){for(j in 1:10){Sigma_sub_inv[i,j]=runif(1,min=-1.2,max=1.2)}}
for(j in 1){for(i in 2:10){Sigma_sub_inv[i,j]=Sigma_sub_inv[j,i]}}
diag(Sigma_sub_inv) <-runif(10,min=3,max=4) 
Sigma_sub <- solve(Sigma_sub_inv)

Sigma <- NULL
for(i in 1:sep){
  temp <- matrix(0,p,p/sep)
  temp[((i-1)*p/sep+1):(i*p/sep),] <- Sigma_sub
  Sigma <- cbind(Sigma,temp)
}

Sigma_inv <- NULL
for(i in 1:sep){
  temp <- matrix(0,p,p/sep)
  temp[((i-1)*p/sep+1):(i*p/sep),] <- Sigma_sub_inv
  Sigma_inv <- cbind(Sigma_inv,temp)
}

#simulate data
X1 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X2 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X3 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X4 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X5 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X6 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X7 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X8 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X9 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X10 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X11 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X12 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X13 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X14 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X15 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X16 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X17 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X18 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X19 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X20 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X <- cbind(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,X18,X19,X20)
#X <- cbind(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10)
#X <- cbind(X1,X2,X3,X4,X5)
#X <- cbind(X1,X2,X3,X4)
X <- cbind(X1,X2)
#X1 <- X[1:200,]
sample_id1 <- sample(1:10,10,replace = FALSE)
sample_id2 <- sample(11:20,0,replace = FALSE)
sample_id <- c(sample_id1,sample_id2)
beta_true1 <- matrix(rep(0, p*2))   #set how many parameters
beta_true1[sample_id] <- rep(1, 10)
#beta_true1[sample_id] <- sample(c(-1,1),10,replace=TRUE)
#combine Sigma_inv for real Sigma_inv match the dimension
Sigma_inv1 <- NULL
for(i in 1:2){
  temp <- matrix(rep(0,p*p*2),nrow=p)
  temp[,((i-1)*p+1):(i*p)] <- Sigma_inv
  Sigma_inv1 <- rbind(Sigma_inv1,temp)
}
Sigma_inv <- Sigma_inv1
beta_true <- Sigma_inv%*%beta_true1    #generate non-zero parameters
#beta_true[which(beta_true1==0),] <- 0
#beta_true <- matrix(runif(p,min=-5,max=5))
p <- exp(X%*%beta_true)/(1+exp(X%*%beta_true))
y <- rbinom(dim(X)[1], 1, p)


#model fitting
#group lasso
glm_mine <- IRLS_IP_gglasso(X=X,y=y,lambda_lasso=0)
#cbind(glm_mine,beta_true)
#sparse group lasso
glm_mine1 <- IRLS_IP_sparsegl(X=X,y=y,lambda_lasso=0)
#cbind(glm_mine1,beta_true)
#glm
glm_r <- glm(y~X,family = binomial) %>% coef
library(glmnet)
#lasso
temp <- glmnet(x=X,y=y, alpha = 1,lambda=0.095, family = binomial(link = "logit"), relax=F,intercept = F)
glm_laaso <- as.matrix(temp[["beta"]])
glm_laaso_bo <- temp[["a0"]]
#length(which(glm_laaso[which(beta_true!=0)] != 0))/length(which(beta_true!=0))
#predict(temp, X, type = "class")

#ridge
temp <- glmnet(x=X,y=y, alpha = 0,lambda=0.05, family = binomial(link = "logit"), relax=F,intercept = F)
glm_ridge <- as.matrix(temp[["beta"]])
glm_ridge_bo <- temp[["a0"]]

#elastic net
temp <- glmnet(x=X,y=y, alpha = 0.5,lambda=0.11, family = binomial(link = "logit"), relax=F,intercept = F)
glm_EN <- as.matrix(temp[["beta"]])
glm_EN_bo <- temp[["a0"]]
#length(which(glm_EN[which(beta_true!=0)] != 0))/length(which(glm_EN!=0))
###
temp <- cbind(glm_mine,glm_mine1,glm_laaso,glm_ridge,glm_EN,glm_r[-1],beta_true) 
colnames(temp) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_ridge","glm_EN","glm","beta_true")
intercept <- cbind(glm_laaso_bo,glm_ridge_bo,glm_EN_bo)
colnames(intercept) <- c("glm_laaso","glm_ridge","glm_EN")
#temp[,c(1,2,3,7)]