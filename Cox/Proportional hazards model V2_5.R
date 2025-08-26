library(survival)
library(ggplot2)
library(dplyr)
library(ranger)
Tol <- 1e-6
ITERS <- 100
INITIAL_BETA <- 0

GetRiskSet <- function(time_of_interest,
                       entry_vector,
                       time_vector,
                       event_vector) {
  return(which((time_of_interest >= entry_vector) & ((time_vector == time_of_interest & event_vector == 1) | (time_vector > time_of_interest))))
  #risk set R(T_j) essentially means that the patient hasn’t failed yet or that their censoring date hasn’t passed yet.
}

GetRiskSet1 <- function(time_of_interest,
                       entry_vector,
                       time_vector,
                       event_vector) {
  return(which((time_vector <= time_of_interest )))
  #risk set R(T_j) essentially means that the patient hasn’t failed yet or that their censoring date hasn’t passed yet.
}

#1. X is a n × p matrix representing p covariates for n patients
#2. x_i is a 1 × p vector representing p covariates for patient j
#3. β is a p × 1 vector representing the coefficients for the p covariates
#Calculating the derivative with respect to β, gradient of the likelihood
CoxGradient <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  p <- ncol(Xs)
  
  gradient <- apply(cbind(Ts, Xs), 1, 
                    function(df){
                      
                      df <- matrix(df, nrow = 1)
                      ts <- df[, 1]
                      xs <- df[, 2:ncol(df)]
                      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event),] %>% 
                        matrix(ncol = ncol(Xs))
                      
                      t1 <- t(X_risk_set) %*% exp(X_risk_set %*% beta)
                      t2 <- sum(exp(X_risk_set %*% beta))
                      
                      return(xs - t1 / t2)
                      
                    }) %>% 
    matrix(nrow = p) %>%
    rowSums()
  return(gradient)
}

#Calculating the Hessian matrix, variance of the likelihood
CoxVariance <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  
  p <- ncol(Xs)
  
  variance <- apply(cbind(Ts, Xs), 1, 
                    function(df){
                      
                      df <- matrix(df, nrow = 1)
                      ts <- df[, 1]
                      xs <- df[, 2:ncol(df)]
                      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event),] %>% 
                        matrix(ncol = ncol(Xs))
                      
                      sum1 <- sum(exp(X_risk_set %*% beta))
                      sum2 <- apply(X_risk_set, 1, 
                                    function(temp){
                                      temp <- matrix(temp, nrow = 1)
                                      temp1 <-  as.numeric(exp(temp%*%beta))* t(temp)%*%temp
                                      return(temp1)
                                    })%>%
                        rowSums() %>%
                        matrix(nrow = p)
                      #sum2 <- rowSums(t(X_risk_set^2) %*% exp(X_risk_set %*% beta))
                      sum3 <- (t(X_risk_set) %*% exp(X_risk_set %*% beta))%*%(t(exp(X_risk_set %*% beta)) %*% X_risk_set)
                      
                      return(- (sum1 * sum2 - sum3) / sum1^2)
                      
                    }) %>%
    rowSums() %>%
    matrix(nrow = p)
  
  return(solve(variance))
  
}

CoxVariance1 <- function(beta,
                        Xs,
                        entry,
                        Ts,
                        event) {
  
  p <- ncol(Xs)
  
  variance <- apply(cbind(Ts, Xs), 1, 
                    function(df){
                      
                      df <- matrix(df, nrow = 1)
                      ts <- df[, 1]
                      xs <- df[, 2:ncol(df)]
                      X_risk_set <- Xs[GetRiskSet(ts, entry, Ts, event),] %>% 
                        matrix(ncol = ncol(Xs))
                      
                      sum1 <- sum(exp(X_risk_set %*% beta))
                      sum2 <- rowSums(t(X_risk_set^2) %*% exp(X_risk_set %*% beta))
                      sum3 <- rowSums(t(X_risk_set) %*% exp(X_risk_set %*% beta))^2
                      
                      return(- (sum1 * sum2 - sum3) / sum1^2)
                      
                    }) %>%
    matrix(nrow = p) %>%
    rowSums()
  
  return(1 / variance)
  
}

Weighted <- function(beta,
                     Xs,
                     entry,
                     Ts,
                     event) {
  p <- ncol(Xs)
  #mu_hat <- beta <- as.matrix(rep(0,nrow(Xs)))
  mu_hat <- apply(cbind(Ts, Xs), 1, 
                  function(df){
                    df <- matrix(df, nrow = 1)
                    ts <- df[, 1]
                    xs <- df[, 2:ncol(df)]
                    t1 <- exp(matrix(xs,nrow = 1)%*%beta)
                    X_risk_set <- Xs[GetRiskSet1(ts, entry, Ts, event),] %>% 
                      matrix(ncol = ncol(Xs))
                    Ts_risk_set <- Ts[GetRiskSet1(ts, entry, Ts, event)]
                    t2 <- apply(cbind(Ts_risk_set, X_risk_set), 1, 
                                function(df1){
                                  df1 <- matrix(df1, nrow = 1)
                                  ts1 <- df[, 1]
                                  xs1 <- df[, 2:ncol(df)] 
                                  t3 <- length(which(Ts_risk_set==ts1))
                                  X_risk_set1 <- Xs[GetRiskSet(ts, entry, Ts, event),] %>% 
                                    matrix(ncol = ncol(Xs))
                                  t4 <- sum(exp(X_risk_set1 %*% beta))
                                  return((t3/t4)*t1)
                                })%>% sum() 
                    
                    return(t2)
                    
                  }) 
  mu_hat <- matrix(mu_hat,ncol=1)
  return(mu_hat)
}

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
    Ni_subset = as.vector(Ni_matrix[,i])
    temp1 = as.numeric(i*Ni_matrix[,i])
    temp1 = temp1[which(temp1!=0)]
    goup_list = c(goup_list,temp1)
    #  X_tilde = X
    #  X_tilde[,which(Ni_subset==0)] <- 0
    temp = X[,which(Ni_subset!=0)]
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
  max((beta0-beta1)^2)<=conv.eps
}

prox_L1 = function(x, threshold) {
  return(sign(x) %*% max(abs(x) - threshold, 0))
}

IRLS_IP_gglasso <- function(x,times,event,entry,lambda_lasso=0,max.iter=100,conv.eps=1e-6){
  store_gradient <- NULL
  store_betas <- matrix(NA, nrow = ITERS, ncol = ncol(x))
  store_variance <- NULL
  beta_prev <- matrix(rep(INITIAL_BETA, ncol(x)), nrow = ncol(x))               #beta_{t-1} (for comparisons)
  
  for(iter in 1:max.iter){
    #iter=1
    if(iter==1){
      store_gradient <- CoxGradient(beta_prev, x, entry, times, event)         #score
      store_variance <- CoxVariance(beta_prev, x, entry, times, event)         #Hessian
      w <- Weighted(beta_prev, x, entry, times, event)
      D <- diag(as.vector(w))                                                  #weighted matrix   
      z <- x%*%beta_prev - x%*%store_variance%*%store_gradient                 #initial response
      z1 <- sqrt(D)%*%z
      X1 <- sqrt(D)%*%x
      #output <- IP_gglasso_calculation(X1,z1,lambda_lasso=0)
      #beta <- output$beta_prediction
      
      lambda_lasso1 = lambda_lasso
      lambda_lasso2 = 0.1
      output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
      beta <- output$beta_prediction
      beta[which(beta<lambda_lasso2&(-lambda_lasso2)<beta)]=0
    } else{
      lambda_lasso1 = 0.0001
      store_gradient <- CoxGradient(beta_prev, x, entry, times, event)
      store_variance <- CoxVariance(beta_prev, x, entry, times, event)
      #store_variance[which(store_variance==Inf)] <- 4000
      w <- Weighted(beta_prev, x, entry, times, event)
      D <- diag(as.vector(w))                                                  #weighted matrix 
      x1 <- x
      x1[,which(beta_prev==0)] <- 0
      z <- x1%*%beta_prev - x1%*%store_variance%*%store_gradient              #initial response
      z1 <- sqrt(D)%*%z
      X1 <- sqrt(D)%*%x1
      if(any(is.nan(z1))=="TRUE"){
        break
      } else {
        output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
        beta <- output$beta_prediction
        beta[which(beta_prev==0)] <- 0
      }
    }
    
    if(did_we_converge(beta_prev,beta,conv.eps)){
      temp <- cbind(abs(beta) - 0,0)
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
