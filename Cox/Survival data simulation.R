
lambdaT = 20 # baseline hazard
lambdaC = 40  # hazard of censoring

library(tidyverse)
library(MASS)
library(matrixcalc)
#p is # of variables
n=100      #sample size
p=40      #parameters size
sep <- 4  # # of x_i in a structure

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
#beta_true1[sample_id] <- rep(1, 10)
beta_true1[sample_id] <- sample(c(-0.5,0.5),10,replace=TRUE)
beta_true <- Sigma_inv%*%beta_true1
# true event time
T = rweibull(n, shape=1, scale=lambdaT*exp(X%*%beta_true)) 
C = rweibull(n, shape=1, scale=lambdaC)   #censoring time
time_raw = pmin(T,C)  #observed time is min of censored and true
event_raw = time_raw==T   # set to 1 if event is observed
table(event_raw)
event_raw <- as.data.frame(event_raw); colnames(event_raw)<-"event"
event_raw[event_raw$event=="TRUE",] <-1
times <- as.matrix(time_raw[which(event_raw$event==1)]); colnames(times)<-"time"
event <- as.matrix(event_raw[which(event_raw$event==1),]) ; colnames(event)<-"event"
x <- as.matrix(X[which(event_raw$event==1),])
entry <- as.matrix(rep(0, dim(times)[1]))

Simulation_data <- as.data.frame(cbind(times,x))
formula <- Surv(time) ~.

