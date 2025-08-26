
lambdaT = 20 # baseline hazard
lambdaC = 35  # hazard of censoring

library(tidyverse)
library(MASS)
library(matrixcalc)
#p is # of variables
n=50     #sample size   50/100/200                                
p=50      #parameters size   200                         Appendix
sep <- 5  # # of x_i in a structure

##############################################T1####################
#Sigma_sub_inv <- matrix(0,p/sep,p/sep)   
#for(i in 1){for(j in 1:10){Sigma_sub_inv[i,j]=runif(1,min=-1.2,max=1.2)}}
#for(j in 1){for(i in 2:10){Sigma_sub_inv[i,j]=Sigma_sub_inv[j,i]}}
#diag(Sigma_sub_inv) <-runif(10,min=3,max=4) 
#Sigma_sub <- solve(Sigma_sub_inv)
##############################################T2####################
Sigma_sub_inv <- matrix(0,p/sep,p/sep)   
for(i in 1:9){Sigma_sub_inv[i,i+1]=runif(1,min=-1.2,max=1.2)}
for(j in 1:9){for(i in 2:10){Sigma_sub_inv[i,j]=Sigma_sub_inv[j,i]}}
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
#simulate data   #change when different sample size,n, and parameter size, p.
X1 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X2 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X3 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X4 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X5 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X6 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X7 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
X8 <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
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
X <- cbind(X1,X2,X3,X4,X5,X6,X7,X8)
#X <- cbind(X1,X2,X3,X4)
#X <- cbind(X1,X2)
#X <- X1

#coefficients
sample_id1 <- sample(1:10,10,replace = FALSE)
sample_id2 <- sample(11:20,0,replace = FALSE)
sample_id <- c(sample_id1,sample_id2)
beta_true1 <- matrix(rep(0, p*8))      #remember to fit the final parameter size
#beta_true1[sample_id] <- rep(1, 10)
beta_true1[sample_id] <- sample(c(-0.75,0.75),10,replace=TRUE)

#beta_true1[sample_id] <- sample(c(-1,1),10,replace=TRUE)
#combine Sigma_inv for real Sigma_inv match the dimension
#remember to fit the final parameter size
Sigma_inv1 <- NULL
for(i in 1:8){
  temp <- matrix(rep(0,p*p*8),nrow=p)
  temp[,((i-1)*p+1):(i*p)] <- Sigma_inv
  Sigma_inv1 <- rbind(Sigma_inv1,temp)
}
Sigma_inv <- Sigma_inv1
beta_true <- Sigma_inv%*%beta_true1 

# true event time
T = rweibull(n, shape=1, scale=lambdaT*exp(-(X%*%beta_true)))               #double check????
#C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
C= quantile(T,probs=0.6)                   #censoring time /study stop time

time_raw = pmin(T,C)  #observed time is min of censored and true
event_raw = time_raw==T   # set to 1 if event is observed
table(event_raw)
event_raw <- as.data.frame(event_raw); colnames(event_raw)<-"event"
event_raw[event_raw$event=="TRUE",] <-1
#times <- as.matrix(time_raw[which(event_raw$event==1)]); colnames(times)<-"time"
data_temp <- cbind(time_raw,event_raw)


###############################################################################
###############################################################################
###############################################################################
#model fitting
#group lasso
survival_mine <- IRLS_IP_gglasso(x=X,times=time_raw,event=event_raw,lambda_lasso=0.07)
#lasso
library(glmnet)
Y_raw <- as.matrix(cbind(time_raw,event_raw)); colnames(Y_raw) <- c("time","status");
temp <- glmnet(x=X, y=Y_raw, family = "cox",lambda=0.14,alpha = 1)
glm_laaso <- as.matrix(temp[["beta"]])
#cbind(beta_true,survival_mine,glm_laaso)
#ridge
temp <- glmnet(x=X, y=Y_raw, family = "cox",lambda=0.15,alpha = 0)
glm_ridge <- as.matrix(temp[["beta"]])
#elastic net
temp <- glmnet(x=X, y=Y_raw, family = "cox",lambda=0.26,alpha = 0.5)
glm_EN <- as.matrix(temp[["beta"]])
###############################################################################
#combine the results
temp <- cbind(survival_mine,glm_laaso,glm_ridge,glm_EN,beta_true) 
colnames(temp) <- c("IRLS_IP_gglasso","cox_lasso","cox_ridge","cox_EN","beta_true")
temp