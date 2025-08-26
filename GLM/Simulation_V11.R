##########################################################################
############################# Simulation #################################
##########################################################################

library(tidyverse)
library(MASS)
library(matrixcalc)
#p is # of variables
n=50     #sample size
p=50      #parameters size
sep <- 5  # # of x_i in a structure

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

X <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X1 <- X[1:200,]
sample_id1 <- sample(1:10,10,replace = FALSE)
sample_id2 <- sample(11:20,0,replace = FALSE)
sample_id <- c(sample_id1,sample_id2)
beta_true1 <- matrix(rep(0, p))
beta_true1[sample_id] <- rep(0.4, 10)
#beta_true1[sample_id] <- sample(c(-1,1),10,replace=TRUE)
beta_true <- Sigma_inv%*%beta_true1
#beta_true[which(beta_true1==0),] <- 0
#beta_true <- matrix(runif(p,min=-5,max=5))
etaTrue <- X%*%beta_true
y <- rpois(dim(X)[1],exp(etaTrue))
max(etaTrue)

#model fitting
#group lasso
try(glm_mine <- IRLS_IP_gglasso(X=X,y=y,lambda_lasso=0))
#cbind(glm_mine,beta_true)
#sparse group lasso
try(glm_mine1 <- IRLS_IP_sparsegl(X=X,y=y,lambda_lasso=0))
#cbind(glm_mine1,beta_true)
#glm
glm_r <- glm(y~X,poisson(link = "log")) %>% coef
library(glmnet)
#lasso
temp <- glmnet(x=X,y=y, alpha = 1,lambda=2, poisson(link = "log"), relax=T,intercept = F)
glm_laaso <- as.matrix(temp[["beta"]])
glm_laaso_bo <- temp[["a0"]]
#length(which(glm_laaso[which(beta_true!=0)] != 0))/length(which(beta_true!=0))
#predict(temp, X, type = "class")

#ridge
temp <- glmnet(x=X,y=y, alpha = 0,lambda=1, poisson(link = "log"), relax=T,intercept = F)
glm_ridge <- as.matrix(temp[["beta"]])
glm_ridge_bo <- temp[["a0"]]

#elastic net
temp <- glmnet(x=X,y=y, alpha = 0.5,lambda=8, poisson(link = "log"), relax=T,intercept = F)
glm_EN <- as.matrix(temp[["beta"]])
glm_EN_bo <- temp[["a0"]]
#length(which(glm_EN[which(beta_true!=0)] != 0))/length(which(glm_EN!=0))
###
temp <- cbind(glm_mine,glm_mine1,glm_laaso,glm_ridge,glm_EN,glm_r[-1],beta_true) 
colnames(temp) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_ridge","glm_EN","glm","beta_true")
intercept <- cbind(glm_laaso_bo,glm_ridge_bo,glm_EN_bo)
colnames(intercept) <- c("glm_laaso","glm_ridge","glm_EN")
#temp[,c(1,2,3,7)]

