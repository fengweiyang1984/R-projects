#Iterative Proximal (IP) Algorithm
#Sparse Regression Incorporating Graphical Structure among Predictors(2016)
#Guan Yu and Yufeng Liu

#Problem is to calculate the proximal beta for the regression fitting
#Matrix X n*p
#Matrix Y n*1

#This is version can successfully running the logit regression
#simulate datasets
setwd("C:/Users/z021w783/Desktop/GRA/Disertation/Results")
#loading function
IP_gglasso_calculation <- function(X,Y,lambda_lasso=0){
#parameters setting
np <- dim(X)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
vnames <- colnames(X)

###### Extract Ni_matrix #########
XtY <- cor(X,Y)
XtX_inv <- round(cor(X),digits = 4)

Ni_matrix = NULL
for (i in 1:nvars) {
  Ni = matrix(rep(1,nvars),nrow = 1)
  Ni[which(XtX_inv[,i]==0)]=0
  Ni[which(as.vector(XtY)==0)]=0
  Ni_matrix = rbind(Ni_matrix,Ni)
}


###### Generate "X~"(X_tilde) by X matrix and Ni_matrix #########
###### Apply gglasso package to generate the "V~"(V_tilde) ######
X_tilde = NULL
goup_list = NULL
for (i in 1:nvars) {
  Ni_subset = as.vector(Ni_matrix[i,])
  temp1 = as.numeric(i*Ni_matrix[i,])
  temp1 = temp1[which(temp1!=0)]
  goup_list = c(goup_list,temp1)
#  X_tilde = X
#  X_tilde[,which(Ni_subset==0)] <- 0
  temp = X[,1:nvars*Ni_subset]
  X_tilde = cbind(X_tilde,temp)
}

fit_ls <- gglasso(x = X_tilde, y = Y, group = goup_list, lambda=lambda_lasso,loss = "ls")
V_Prediction <- fit_ls$beta
V_slide_matrix = matrix(rep(0,nvars*nvars),nrow = nvars)
for (i in 1:nvars) {
  Ni_subset = as.vector(Ni_matrix[,i])
  V_slide_matrix[which(Ni_subset!=0),i]=V_Prediction[which(goup_list==i)]
}
beta_prediction = as.matrix(apply(V_slide_matrix,1,sum))
X_selected = X[beta_prediction!=0]
Y_prediction = X%*%beta_prediction
output = list(X_selected,beta_prediction,Y_prediction,V_slide_matrix)
names(output) = c("X_selected","beta_prediction","Y_prediction","V_slide_matrix")
return(output)
}

did_we_converge <- function(beta0,beta1,conv.eps){
  sum(abs(beta0-beta1))<=conv.eps
}

prox_L1 = function(x, threshold) {
  return(sign(x) %*% max(abs(x) - threshold, 0))
}

IRLS_IP_gglasso <- function(X,y,lambda_lasso=0,max.iter=100,conv.eps=1e-10){
  library(gglasso)
  beta <- as.matrix(rep(0,ncol(X))) #initialize beta
  beta_prev <- beta               #beta_{t-1} (for comparisons)
  
  for(iter in 1:max.iter){
    #iter=1
      if(iter==1){
        for (i in 1:length(y)){
          if (y[i]==0) y[i] =1e-16}
        mu <- as.matrix(y)
        eta <- log(mu)
        z <- eta + (y-mu)/mu                  #initial response
        w <- c(mu)  
        
        D <- diag(w)
        z1 <- sqrt(D)%*%z
        X1 <- sqrt(D)%*%X
        lambda_lasso1 = 0.1
        output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
        beta <- output$beta_prediction

        #p[which(X%*%beta>40)] <- 0.999999
        #p[which(X%*%beta< (-40))] <- 10e-8
        beta[which(beta<0.04&(-0.04)<beta)]=0
        beta[which(beta>lambda_lasso1)]=beta[which(beta>lambda_lasso1)]+lambda_lasso1
        beta[which(beta<(-lambda_lasso1))]=beta[which(beta<(-lambda_lasso1))]-lambda_lasso1
      } else{
        lambda_lasso1 = 0.0007
        X1 <- X[,which(beta_prev!=0)]
        beta1 <- matrix(beta[which(beta!=0)])
        eta<-X1%*%beta1
        z=eta+(y-exp(eta))/exp(eta)
        #p[which(X1%*%beta1>40)] <- 0.999999
        #p[which(X1%*%beta1< (-40))] <- 10e-8
        V<-c(exp(eta))
        
        D <- diag(V)
        z1 <- sqrt(D)%*%z
        X1 <- sqrt(D)%*%X1
        output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
        beta_temp <- output$beta_prediction
        #beta[which(beta>0)]=beta[which(beta>0)]+0.03
        #beta[which(beta<0)]=beta[which(beta<0)]-0.03
        beta[which(beta_prev==0)] <- 0
        beta[which(beta_prev!=0)] <- beta_temp
        #beta[which(beta>0)]=beta[which(beta>0)]+lambda_lasso1
        #beta[which(beta<0)]=beta[which(beta<0)]-lambda_lasso1  
      }

    if(did_we_converge(beta_prev,beta,conv.eps)){
      #temp <- cbind(abs(beta) - lambda_lasso,0)
      #temp <- apply(temp,1,max)
      #beta <- sign(beta) *temp
      break
    } else {
      beta_prev <- beta
    }
    print(cbind(beta,beta_prev,beta_true))
    print(iter)
  }
  return(beta)
  
}
#glm_mine <- IRLS_IP_gglasso(X=X,y=y)

##########################################################################
############################# Simulation #################################
##########################################################################
library(tidyverse)
library(MASS)
#p is # of variables
n=50
p=50
sep <- 5
Sigma_sub <- matrix(NA,p/sep,p/sep)   
for(i in 1:p/sep){for(j in 1:p/sep){Sigma_sub[i,j]=0.6^(abs(i-j))}}
Sigma <- NULL
for(i in 1:sep){
  temp <- matrix(0,p,p/sep)
  temp[((i-1)*p/sep+1):(i*p/sep),] <- Sigma_sub
  Sigma <- cbind(Sigma,temp)
}
#simulate data
X <- mvrnorm(n = n, rep(0, p), Sigma,empirical = TRUE)
#X1 <- X[1:200,]
sample_id1 <- sample(1:10,5,replace = FALSE)
sample_id2 <- sample(11:20,5,replace = FALSE)
sample_id <- c(sample_id1,sample_id2)
beta_true <- matrix(rep(0, p))
#beta_true[sample_id] <- rnorm(10,0,1)
beta_true[sample_id] <- sample(union(runif(20,min=-2,max=-1),runif(20,min=1,max=2)),10,replace = FALSE)
#beta_true <- matrix(runif(p,min=-5,max=5))
etaTrue <- X%*%beta_true
y <- rpois(dim(X)[1],exp(etaTrue))
max(etaTrue)

#model fitting
glm_mine <- IRLS_IP_gglasso(X=X,y=y,lambda_lasso=0)
glm_r <- glm(y~X,family = poisson) %>% coef
temp <- cbind(glm_mine,glm_r[-1],beta_true)
colnames(temp) <- c("IRLS_IP_gglasso","glm","beta_true")
temp
write.csv(temp,"temp.csv")
