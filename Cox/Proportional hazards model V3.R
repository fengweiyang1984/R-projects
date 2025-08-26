library(survival)
library(ggplot2)
library(dplyr)
library(ranger)

fid <- function (x, index) 
{
  idup = duplicated(x)
  if (!any(idup)) 
    list(index_first = index, index_ties = NULL)
  else {
    ndup = !idup
    xu = x[ndup]
    index_first = index[ndup]
    ities = match(x, xu)
    index_ties = split(index, ities)
    nties = sapply(index_ties, length)
    list(index_first = index_first, index_ties = index_ties[nties > 
                                                              1])
  }
}

coxgrad <- function (eta, y, w, std.weights = FALSE, diag.hessian = TRUE) 
{
  nobs <- nrow(y)
  if (missing(w)) 
    w = rep(1, length(eta))
  if (std.weights) 
    w = w/sum(w)

  if ("strata" %in% names(attributes(y))) {
    strata <- attr(y, "strata")
  }
  else {
    strata <- rep(1, nobs)
  }
  if (length(strata) != nobs) 
    stop("length of strata != nobs")
  if (length(unique(strata)) == 1) {
    time <- y[, "times"]
    d <- y[, "event"]
    eta <- scale(eta, TRUE, FALSE)
    if ("stop_time" %in% names(attributes(y))) {
      o <- attr(y, "stop_time")
    }
    else {
      o <- order(time, d, decreasing = c(FALSE, TRUE))
    }
    exp_eta <- exp(eta)[o]
    time <- time[o]
    d <- d[o]
    w <- w[o]
    rskden <- rev(cumsum(rev(exp_eta * w)))
    dups <- fid(time[d == 1], seq(length(d))[d == 1])
    dd <- d
    ww <- w
    if (!is.null(ties <- dups$index_ties)) {
      dd[unlist(ties)] = 0
      dd[dups$index_first] = 1
      wsum = sapply(ties, function(i, w) sum(w[i]), ww)
      tie1 = sapply(ties, function(i) i[1])
      ww[tie1] = wsum
    }
    rskcount = cumsum(dd)
    rskdeninv = cumsum((ww/rskden)[dd == 1])
    rskdeninv = c(0, rskdeninv)
    grad <- w * (d - exp_eta * rskdeninv[rskcount + 1])
    grad[o] <- grad
    if (diag.hessian) {
      rskdeninv2 <- cumsum((ww/(rskden^2))[dd == 1])
      rskdeninv2 <- c(0, rskdeninv2)
      w_exp_eta <- w * exp_eta
      diag_hessian <- w_exp_eta^2 * rskdeninv2[rskcount + 
                                                 1] - w_exp_eta * rskdeninv[rskcount + 1]
      diag_hessian[o] <- diag_hessian
      attr(grad, "diag_hessian") <- diag_hessian
    }
    return(grad)
  }
  else {
    overall_grad <- rep(NA, nobs)
    if (diag.hessian) 
      overall_diag_hessian <- rep(NA, nobs)
    for (i in unique(strata)) {
      ii <- which(strata == i)
      strata_res <- coxgrad2(eta[ii], y[ii, , drop = FALSE], 
                             w[ii], std.weights = FALSE, diag.hessian = diag.hessian)
      overall_grad[ii] <- strata_res
      if (diag.hessian) {
        overall_diag_hessian[ii] <- attr(strata_res, 
                                         "diag_hessian")
      }
    }
    if (diag.hessian) {
      attr(overall_grad, "diag_hessian") <- overall_diag_hessian
    }
    return(overall_grad)
  }
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
  sum((beta0-beta1)^2)<=conv.eps
}

prox_L1 = function(x, threshold) {
  return(sign(x) %*% max(abs(x) - threshold, 0))
}

IRLS_IP_gglasso <- function(x,times,event,lambda_lasso=0,max.iter=100,conv.eps=1e-6,INITIAL_BETA=0){
  store_gradient <- NULL
  store_betas <- matrix(NA, nrow = max.iter, ncol = ncol(x))
  store_variance <- NULL
  beta_prev <- matrix(rep(INITIAL_BETA, ncol(x)), nrow = ncol(x))               #beta_{t-1} (for comparisons)
  offset = rep(0, nrow(x))
  
  for(iter in 1:max.iter){
    #iter=1
    if(iter==1){
      eta <- x%*%beta_prev
      coxgrad_results <- coxgrad(eta=eta, y=cbind(times,event), diag.hessian = TRUE)
      w <- -attributes(coxgrad_results)$diag_hessian               #weighted matrix
      z <- (eta - offset) - ifelse(w != 0, -coxgrad_results / w, 0)    #response
      D <- diag(as.vector(w))                                                    
      z1 <- sqrt(D)%*%z
      X1 <- sqrt(D)%*%x
      #output <- IP_gglasso_calculation(X1,z1,lambda_lasso=0)
      #beta <- output$beta_prediction
      
      lambda_lasso1 = lambda_lasso
      lambda_lasso2 = 0.001
      output <- IP_gglasso_calculation(X1,z1,lambda_lasso=lambda_lasso1)
      beta <- output$beta_prediction
      beta[which(beta<lambda_lasso2&(-lambda_lasso2)<beta)]=0
    } else{
      lambda_lasso1 =lambda_lasso/4
      eta <- x%*%beta_prev
      coxgrad_results <- coxgrad(eta=eta, y=cbind(times,event), diag.hessian = TRUE)
      w <- -attributes(coxgrad_results)$diag_hessian   #weighted matrix
      z <- (eta - offset) - ifelse(w != 0, -coxgrad_results / w, 0)  #response                           
      D <- diag(as.vector(w))                                                    
      z1 <- sqrt(D)%*%z
      x1 <- x
      x1[,which(beta_prev==0)] <- 0
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
