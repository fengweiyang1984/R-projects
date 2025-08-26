#Iterative Proximal (IP) Algorithm
#Sparse Regression Incorporating Graphical Structure among Predictors(2016)
#Guan Yu and Yufeng Liu

#Problem is to calculate the proximal beta for the regression fitting
#Matrix X n*p
#Matrix Y n*1

#This is version can successfully running the logit regression
#simulate datasets
#setwd("C:/Users/z021w783/Desktop/GRA/Disertation/Results")
#loading function
IP_gglasso_calculation <- function(X,Y,lambda_lasso=0){
  library(gglasso)
  #parameters setting
np <- dim(X)
nobs <- as.integer(np[1])
nvars <- as.integer(np[2])
vnames <- colnames(X)

###### Extract Ni_matrix #########
XtY <- rep(1,nobs)
XtX_inv <- Sigma_inv

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
  temp = X[,which(Ni_subset!=0)]
  X_tilde = cbind(X_tilde,temp)
}

fit_ls <- gglasso(x = X_tilde, y = Y, group = goup_list, lambda=lambda_lasso,loss = "ls",
                  maxit = 10000,eps = 1e-04)
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
  sum((beta0-beta1)^2)<=conv.eps
}

prox_L1 = function(x, threshold) {
  return(sign(x) %*% max(abs(x) - threshold, 0))
}

IRLS_SUR_gglasso <- function(X,y,lambda_lasso=0,max.iter=100,conv.eps=1e-6){
  beta <- as.matrix(rep(0,ncol(X))) #initialize beta
  beta_prev <- beta               #beta_{t-1} (for comparisons)
  
  for(iter in 1:max.iter){
    #iter=1
      if(iter==1){
        temp <- y
        for (i in 1:length(temp)){
          if (temp[i]==0) temp[i] =1e-16}
        mu <- as.matrix(temp)
        eta <- log(mu)
        z <- eta + (y-mu)/mu                  #initial response
        w <- c(mu)  
        D <- diag(w)
        z1 <- sqrt(D)%*%z
        X1 <- sqrt(D)%*%X
        lambda_lasso1 = 3
        lambda_lasso2 = 0.15
        output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
        beta <- output$beta_prediction
        beta[which(beta<lambda_lasso2&(-lambda_lasso2)<beta)]=0
      } else{
        lambda_lasso1 = 0.5
        X1 <- X
        X1[,which(beta_prev==0)] <- 0
        #X1 <- as.matrix(X1)
        beta1 <- beta
        eta<-X1%*%beta1
        z=eta+(y-exp(eta))/exp(eta)
        #p[which(X1%*%beta1>40)] <- 0.999999
        #p[which(X1%*%beta1< (-40))] <- 10e-8
        V<-c(exp(eta))
        z1 <- X1%*%beta1+(y-p)/V
        D <- diag(V)
        z1 <- sqrt(D)%*%z
        X1 <- sqrt(D)%*%X1
        if(any(is.nan(z1))=="TRUE"){
          break
        } else {
          output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
          beta_temp <- output$beta_prediction
          beta <- beta_temp
          beta[which(beta_prev==0)] <-0
        }
      }

    if(did_we_converge(beta_prev,beta,conv.eps)){
      temp <- cbind(abs(beta) - lambda_lasso,0)
      temp <- apply(temp,1,max)
      beta <- sign(beta) *temp
      break
    } else {
      beta_prev <- beta
    }
    #print(cbind(beta,beta_prev,beta_true))
    print(iter)
  }
  return(beta)
  
}
#glm_mine <- IRLS_IP_gglasso(X=X,y=y)

