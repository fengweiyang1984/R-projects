##########################################################################
############################# Simulation Resutl Summary ##################
##########################################################################
setwd("C:/Users/z021w783/Desktop/GRA/Disertation/Results")
#creat tables
temp_result <- matrix(rep(0,11), nrow =1, byrow=TRUE)
colnames(temp_result) = c('false_p','false_n','sensitivity','specificity','MSE','SSE','MAE','false_p','false_n','sensitivity','specificity')
summary_result_1<- NULL
summary_result_2<- NULL
summary_result_3<- NULL
summary_result_4<- NULL
summary_result_5<- NULL
summary_result_6<- NULL
T <- 0

#loading GLM IP methods of gglasso and sparsegl
repeat{

  for (j in 1:5) {
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/GLM Iterative Proximal (IP) Algorithm Sparse Regression for logistic V_gglasso2.R")
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/GLM Iterative Proximal (IP) Algorithm Sparse Regression for logistic V_sparsegl2.R")
#Implement simulation
source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/Simulation.R")
#source("C:/Users/z021w783/Desktop/GRA/Disertation/Code/GLM Paper final library/SimulationV2.R")
temp
  
#gglasso_result <- 1
############################################################
  T1 <- as.vector(temp[,1])
  T2 <- as.vector(temp[,7])
  temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
  temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
  temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T1 != 0)) #sensitivity for variable selection
  temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T1 == 0)) #specificity for variable selection
  temp_result[1,5] <- mean((T1-T2)^2)  #MSE
  temp_result[1,6] <- sum((T1-T2)^2)  #SSE
  temp_result[1,7] <- mean(abs(T1-T2))  #MAE
  p_prediction <- exp(X%*%T1)/(1+exp(X%*%T1))
  y_prediction <- rbinom(dim(X)[1], 1, p_prediction)
  temp_result[1,8] <- length(which(y_prediction[which(y==0)] != 0)) #false positive for outcome prediction
  temp_result[1,9] <- length(which(y_prediction[which(y!=0)] == 0)) #false negative for outcome prediction
  temp_result[1,10] <- length(which(y_prediction[which(y!=0)] != 0))/length(which(y_prediction != 0)) #sensitivity for outcome prediction
  temp_result[1,11] <- length(which(y_prediction[which(y==0)] == 0))/length(which(y_prediction == 0)) #specificity for outcome prediction
  summary_result_1 <- rbind(temp_result,summary_result_1)
  
#sparsegl_result <- 2
  ############################################################
  T1 <- as.vector(temp[,2])
  T2 <- as.vector(temp[,7])
  temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
  temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
  temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T1 != 0)) #sensitivity for variable selection
  temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T1 == 0)) #specificity for variable selection
  temp_result[1,5] <- mean((T1-T2)^2)  #MSE
  temp_result[1,6] <- sum((T1-T2)^2)  #SSE
  temp_result[1,7] <- mean(abs(T1-T2))  #MAE
  p_prediction <- exp(X%*%T1)/(1+exp(X%*%T1))
  y_prediction <- rbinom(dim(X)[1], 1, p_prediction)
  temp_result[1,8] <- length(which(y_prediction[which(y==0)] != 0)) #false positive for outcome prediction
  temp_result[1,9] <- length(which(y_prediction[which(y!=0)] == 0)) #false negative for outcome prediction
  temp_result[1,10] <- length(which(y_prediction[which(y!=0)] != 0))/length(which(y_prediction != 0)) #sensitivity for outcome prediction
  temp_result[1,11] <- length(which(y_prediction[which(y==0)] == 0))/length(which(y_prediction == 0)) #specificity for outcome prediction
  summary_result_2 <- rbind(temp_result,summary_result_2)  
  
#lasso_result <- 3
  ############################################################
  T1 <- as.vector(temp[,3])
  T2 <- as.vector(temp[,7])
  temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
  temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
  temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T1 != 0)) #sensitivity for variable selection
  temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T1 == 0)) #specificity for variable selection
  temp_result[1,5] <- mean((T1-T2)^2)  #MSE
  temp_result[1,6] <- sum((T1-T2)^2)  #SSE
  temp_result[1,7] <- mean(abs(T1-T2))  #MAE
  p_prediction <- exp(X%*%T1)/(1+exp(X%*%T1))
  y_prediction <- rbinom(dim(X)[1], 1, p_prediction)
  temp_result[1,8] <- length(which(y_prediction[which(y==0)] != 0)) #false positive for outcome prediction
  temp_result[1,9] <- length(which(y_prediction[which(y!=0)] == 0)) #false negative for outcome prediction
  temp_result[1,10] <- length(which(y_prediction[which(y!=0)] != 0))/length(which(y_prediction != 0)) #sensitivity for outcome prediction
  temp_result[1,11] <- length(which(y_prediction[which(y==0)] == 0))/length(which(y_prediction == 0)) #specificity for outcome prediction
  summary_result_3 <- rbind(temp_result,summary_result_3)  
  
#ridge_result <- 4
  ############################################################
  T1 <- as.vector(temp[,4])
  T2 <- as.vector(temp[,7])
  temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
  temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
  temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T1 != 0)) #sensitivity for variable selection
  temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T1 == 0)) #specificity for variable selection
  temp_result[1,5] <- mean((T1-T2)^2)  #MSE
  temp_result[1,6] <- sum((T1-T2)^2)  #SSE
  temp_result[1,7] <- mean(abs(T1-T2))  #MAE
  p_prediction <- exp(X%*%T1)/(1+exp(X%*%T1))
  y_prediction <- rbinom(dim(X)[1], 1, p_prediction)
  temp_result[1,8] <- length(which(y_prediction[which(y==0)] != 0)) #false positive for outcome prediction
  temp_result[1,9] <- length(which(y_prediction[which(y!=0)] == 0)) #false negative for outcome prediction
  temp_result[1,10] <- length(which(y_prediction[which(y!=0)] != 0))/length(which(y_prediction != 0)) #sensitivity for outcome prediction
  temp_result[1,11] <- length(which(y_prediction[which(y==0)] == 0))/length(which(y_prediction == 0)) #specificity for outcome prediction
  summary_result_4 <- rbind(temp_result,summary_result_4)
  
#EN_result <- 5
  ############################################################
  T1 <- as.vector(temp[,5])
  T2 <- as.vector(temp[,7])
  temp_result[1,1] <- length(which(T1[which(T2==0)] != 0)) #false positive for variable selection
  temp_result[1,2] <- length(which(T1[which(T2!=0)] == 0)) #false negative for variable selection
  temp_result[1,3] <- length(which(T1[which(T2!=0)] != 0))/length(which(T1 != 0)) #sensitivity for variable selection
  temp_result[1,4] <- length(which(T1[which(T2==0)] == 0))/length(which(T1 == 0)) #specificity for variable selection
  temp_result[1,5] <- mean((T1-T2)^2)  #MSE
  temp_result[1,6] <- sum((T1-T2)^2)  #SSE
  temp_result[1,7] <- mean(abs(T1-T2))  #MAE
  p_prediction <- exp(X%*%T1)/(1+exp(X%*%T1))
  y_prediction <- rbinom(dim(X)[1], 1, p_prediction)
  temp_result[1,8] <- length(which(y_prediction[which(y==0)] != 0)) #false positive for outcome prediction
  temp_result[1,9] <- length(which(y_prediction[which(y!=0)] == 0)) #false negative for outcome prediction
  temp_result[1,10] <- length(which(y_prediction[which(y!=0)] != 0))/length(which(y_prediction != 0)) #sensitivity for outcome prediction
  temp_result[1,11] <- length(which(y_prediction[which(y==0)] == 0))/length(which(y_prediction == 0)) #specificity for outcome prediction
  summary_result_5 <- rbind(temp_result,summary_result_5)
  
#EN_result <- 6
  ############################################################

}
  T = T+1
  if (T == 15){
    break
  }
}

write.csv(summary_result_1,file=paste('summary_result',1,'.csv',seq=""))
write.csv(summary_result_2,file=paste('summary_result',2,'.csv',seq=""))
write.csv(summary_result_3,file=paste('summary_result',3,'.csv',seq=""))
write.csv(summary_result_4,file=paste('summary_result',4,'.csv',seq=""))
write.csv(summary_result_5,file=paste('summary_result',5,'.csv',seq=""))
write.csv(summary_result_6,file=paste('summary_result',6,'.csv',seq=""))

colnames(temp)
