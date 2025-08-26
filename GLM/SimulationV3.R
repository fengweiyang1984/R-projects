##########################################################################
############################# Simulation #################################
##########################################################################

library(tidyverse)
library(MASS)
#p is # of variables
n=200
p=200
sep <- 40
Sigma_sub <- matrix(NA,p/sep,p/sep)   
for(i in 1:p/sep){for(j in 1:p/sep){Sigma_sub[i,j]=0.6^(abs(i-j))}}
Sigma <- NULL
for(i in 1:sep){
  temp <- matrix(0,p,p/sep)
  temp[((i-1)*p/sep+1):(i*p/sep),] <- Sigma_sub
  Sigma <- cbind(Sigma,temp)
}

n1=200
p1=100
sep <- 20
Sigma_sub1 <- matrix(NA,p1/sep,p1/sep)   
for(i in 1:p1/sep){for(j in 1:p1/sep){Sigma_sub1[i,j]=0.6^(abs(i-j))}}
Sigma1 <- NULL
for(i in 1:sep){
  temp <- matrix(0,p1,p1/sep)
  temp[((i-1)*p1/sep+1):(i*p1/sep),] <- Sigma_sub1
  Sigma1 <- cbind(Sigma1,temp)
}
#simulate data
X1 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X2 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X3 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X4 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X5 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X3 <- mvrnorm(n = n1, rep(0, p1), Sigma1,empirical = TRUE)
X <- cbind(X1,X2,X3)
#X1 <- X[1:200,]
sample_id1 <- sample(1:10,8,replace = FALSE)
sample_id2 <- sample(11:20,9,replace = FALSE)
sample_id <- c(sample_id1,sample_id2)
beta_true <- matrix(rep(0, p*2+p1))
beta_true[sample_id] <- sample(union(runif(20,min=-5,max=-3),runif(20,min=3,max=5)),17,replace = FALSE)
#beta_true <- matrix(runif(p,min=-5,max=5))
p <- exp(X%*%beta_true)/(1+exp(X%*%beta_true))
y <- rbinom(dim(X)[1], 1, p)


#model fitting
#group lasso
glm_mine <- IRLS_IP_gglasso(X=X,y=y,lambda_lasso=0)
cbind(glm_mine,beta_true)
#sparse group lasso
glm_mine1 <- IRLS_IP_sparsegl(X=X,y=y,lambda_lasso=0)
cbind(glm_mine1,beta_true)
#glm
glm_r <- glm(y~X,family = binomial) %>% coef
library(glmnet)
#lasso
temp <- glmnet(x=X,y=y, alpha = 1,lambda=0.03, family = binomial(link = "logit"), relax=TRUE)
glm_laaso <- as.matrix(temp[["beta"]])
length(which(glm_laaso[which(beta_true!=0)] != 0))/length(which(glm_laaso!=0))
#ridge
temp <- glmnet(x=X,y=y, alpha = 0,lambda=0.03, family = binomial(link = "logit"), relax=TRUE)
glm_ridge <- as.matrix(temp[["beta"]])
#elestic net
temp <- glmnet(x=X,y=y, alpha = 0.5,lambda=0.03, family = binomial(link = "logit"), relax=TRUE)
glm_EN <- as.matrix(temp[["beta"]])
###
temp <- cbind(glm_mine,glm_mine1,glm_laaso,glm_ridge,glm_EN,glm_r[-1],beta_true) 
colnames(temp) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_ridge","glm_EN","glm","beta_true")

