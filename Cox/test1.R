time <- c(1, 3, 5, 6, 2, 7, 9, 11)
status <- c(1, 0, 1, 1, 1, 0, 1, 1)
sex <- rpois(8,10)
age <- c(57, 52, 48, 42, 39, 31, 26, 22)

df <- data.frame(time, status, sex, age)
library(survival)
fit_cox <- coxph(Surv(time, status) ~ sex + age, data=df, method = "breslow")

breslow_est_adj <- function(time, status, X, B){
  data <- data.frame(time,status,X)
  ID <- 1:dim(X)[1]
  data <- cbind(ID,data)
  colnames(data)[1:3] <- c("ID","time","status")
  data <- data[order(data$time), ]
  data1 <- data[!duplicated(data$time), ]
  t   <- unique(data1$time)
  k    <- length(t)
  h    <- rep(0,k)
    
  for(i in 1:k) {
     LP_sample <- sum(colMeans(X) * B) 
     LP_indiv <- c(as.matrix(data[,-1:-3])%*%as.matrix(B)) 
     lp_centered <- (LP_indiv - LP_sample)[data$time>=t[i]]
     risk <- exp(lp_centered)
     h[i] <- sum(data$status[data$time==t[i]]) / sum(risk)
     }
   
  H0 <- cumsum(h)
  H0 <- cbind(as.data.frame(H0),data1)
  H01 <- as.data.frame(data[,1:3])
  library(dplyr) 
  H01 <- full_join(H01, H0[,1:4], by=c("time"))
  H01 <- H01[,1:4]; colnames(H01) <- c("ID","time","status","H0")
  return(H01)
  }


Survival_Prob_est<- function(time2, status2, X2, B2){
  H0 <- breslow_est_adj(time=time2, status=status2, X=X2, B=B2)
  data <- data.frame(time2,status2,X2)
  ID <- 1:dim(X)[1]
  data <- cbind(ID,data)
  #LP <- predict(fit_cox, type="lp")  
  LP_indiv <- c(as.matrix(data[,-1:-3])%*%as.matrix(B2)) 
  LP_sample <- sum(colMeans(X2) * B2) 
  LP <- LP_indiv - LP_sample
  res <- exp(-H0[match(data$ID, H0$ID),4]*exp(LP))
  return(res)
}


temp <- Survival_Prob_est(time2=df$time, status2=df$status, X2=df[,-1:-2], B2=fit_cox$coef)
temp

fit_cox <- coxph(Surv(time, status) ~ sex + age, data=df, method = "breslow")
predict(fit_cox, type="survival") 

H0 <- basehaz(fit_cox,centered = TRUE) 
H0
LP <- predict(fit_cox, type="survival") 
LP
exp(-H0[, 1]*exp(LP))
