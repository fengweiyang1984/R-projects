##########################################################################
############################# Simulation Resutl Summary ##################
##########################################################################

#creat tables
setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Results")
#creat tables
temp_result <- matrix(rep(0,9), nrow =1, byrow=TRUE)
colnames(temp_result) = c('false_p','false_n','sensitivity','specificity','MSE','SSE','MAE','Brier','C-index')
summary_result_1<- NULL
summary_result_2<- NULL
summary_result_3<- NULL
summary_result_4<- NULL
RT <- 0

#loading GLM IP methods of gglasso and sparsegl
repeat{
  
  for (j in 1:10) {
    #Implement simulation
    source("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Code/Liabary/Proportional hazards model V3_1.R")    
    source("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Code/Liabary/Survival data simulation V1_3.R")
    library(survival)
    library(SurvMetrics)
    temp_event <- sapply(event_raw,as.numeric)
    
    #gglasso_result <- 1
    ############################################################
    T1 <- as.vector(temp[,1])
    T2 <- as.vector(temp[,5])
    temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
    temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
    temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T2 != 0)) #sensitivity for variable selection
    temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T2 == 0)) #specificity for variable selection
    temp_result[1,5] <- mean((T1-T2)^2)  #MSE
    temp_result[1,6] <- sum((T1-T2)^2)  #SSE
    temp_result[1,7] <- mean(abs(T1-T2))  #MAE
    temp_SurvProb <- Survival_Prob_est(time2=time_raw, status2=event_raw, X2=X, B2=T1)
    temp_result[1,8] <- Brier(Surv(time=time_raw, event=as.numeric(temp_event)), temp_SurvProb)
    cindex_temp <- cindex_calcu(time2=time_raw, status2=temp_event, X2=X, B2=T1)[[3]]
    temp_result[1,9] <- cindex_temp
    summary_result_1 <- rbind(temp_result,summary_result_1)
    
    #lasso_result <- 2
    ############################################################
    T1 <- as.vector(temp[,2])
    T2 <- as.vector(temp[,5])
    temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
    temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
    temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T2 != 0)) #sensitivity for variable selection
    temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T2 == 0)) #specificity for variable selection
    temp_result[1,5] <- mean((T1-T2)^2)  #MSE
    temp_result[1,6] <- sum((T1-T2)^2)  #SSE
    temp_result[1,7] <- mean(abs(T1-T2))  #MAE
    temp_SurvProb <- Survival_Prob_est(time2=time_raw, status2=event_raw, X2=X, B2=T1)
    temp_result[1,8] <- Brier(Surv(time=time_raw, event=as.numeric(temp_event)), temp_SurvProb)
    cindex_temp <- cindex_calcu(time2=time_raw, status2=temp_event, X2=X, B2=T1)[[3]]
    temp_result[1,9] <- cindex_temp
    summary_result_2 <- rbind(temp_result,summary_result_2)  

    #ridge_result <- 3
    ############################################################
    T1 <- as.vector(temp[,3])
    T2 <- as.vector(temp[,5])
    temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
    temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
    temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T2 != 0)) #sensitivity for variable selection
    temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T2 == 0)) #specificity for variable selection
    temp_result[1,5] <- mean((T1-T2)^2)  #MSE
    temp_result[1,6] <- sum((T1-T2)^2)  #SSE
    temp_result[1,7] <- mean(abs(T1-T2))  #MAE
    temp_SurvProb <- Survival_Prob_est(time2=time_raw, status2=event_raw, X2=X, B2=T1)
    temp_result[1,8] <- Brier(Surv(time=time_raw, event=as.numeric(temp_event)), temp_SurvProb)
    cindex_temp <- cindex_calcu(time2=time_raw, status2=temp_event, X2=X, B2=T1)[[3]]
    temp_result[1,9] <- cindex_temp
    summary_result_3 <- rbind(temp_result,summary_result_3)  
    
    #EN_result <- 4
    ############################################################
    T1 <- as.vector(temp[,4])
    T2 <- as.vector(temp[,5])
    temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
    temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
    temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T2 != 0)) #sensitivity for variable selection
    temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T2 == 0)) #specificity for variable selection
    temp_result[1,5] <- mean((T1-T2)^2)  #MSE
    temp_result[1,6] <- sum((T1-T2)^2)  #SSE
    temp_result[1,7] <- mean(abs(T1-T2))  #MAE
    temp_SurvProb <- Survival_Prob_est(time2=time_raw, status2=event_raw, X2=X, B2=T1)
    temp_result[1,8] <- Brier(Surv(time=time_raw, event=as.numeric(temp_event)), temp_SurvProb)
    cindex_temp <- cindex_calcu(time2=time_raw, status2=temp_event, X2=X, B2=T1)[[3]]
    temp_result[1,9] <- cindex_temp
    summary_result_4 <- rbind(temp_result,summary_result_4)  

  }
  RT <- RT+1
  if (RT == 15){
    break
  }
}

write.csv(summary_result_1,file=paste('summary_result',1,'.csv',seq=""))
write.csv(summary_result_2,file=paste('summary_result',2,'.csv',seq=""))
write.csv(summary_result_3,file=paste('summary_result',3,'.csv',seq=""))
write.csv(summary_result_4,file=paste('summary_result',4,'.csv',seq=""))

temp_result<-NULL

temp <- apply(summary_result_1,2,mean) 
temp_result <- rbind(temp_result,temp)
temp <- apply(summary_result_2,2,mean) 
temp_result <- rbind(temp_result,temp)
temp <- apply(summary_result_3,2,mean) 
temp_result <- rbind(temp_result,temp)
temp <- apply(summary_result_4,2,mean) 
temp_result <- rbind(temp_result,temp)
rownames(temp_result) <- c("IRLS_gglasso","lasso","ridge","Elastic Net")
temp_result

