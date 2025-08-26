###################################################
############# Analysis Codes to Predict Test Data Outcomes 
############# 09/06/2020
############# Yanming Li
#############

setwd("C:\\Users\\fengw\\OneDrive\\×ÀÃæ\\GRA-20200830T211630Z-001\\KUCC Research Symposium WEEK 2020")

dat_Train <- as.matrix(readRDS("data_train_10X.rds"))
dim(dat_Train)
#[1] 33538  6292
table(dat_Train[,1])

dat_Test <- as.matrix(readRDS("data_test_10X.rds"))
dim(dat_Test)

y_Train <- as.matrix(readRDS("truth_train_10X.rds"))
y_Train <- as.vector(y_Train[,1])

table(y_Train)

write.csv(dat_Train, file="data_train_10X.csv", row.names=TRUE)
write.csv(dat_Test, file="data_test_10X.csv", row.names=TRUE)
write.csv(y_Train, file="y_Train_10X.csv", row.names=FALSE)

#####################################
#####################################

setwd("C:\\Users\\yli8\\Desktop\\SingleCellClassification\\Kong1")

dat_Train <- as.matrix(readRDS("data_train_Kong.rds"))
dim(dat_Train)
#[1] 33538  1098
table(dat_Train[,1])
dat_Test <- as.matrix(readRDS("data_test_Kong.rds"))
dim(dat_Test)

y_Train <- as.matrix(readRDS("truth_train_Kong.rds"))
y_Train <- as.vector(y_Train[,1])

table(y_Train)

write.csv(dat_Train, file="data_train_Kong.csv", row.names=TRUE)
write.csv(dat_Test, file="data_test_Kong.csv", row.names=TRUE)
write.csv(y_Train, file="y_Train_Kong.csv", row.names=FALSE)

#####################################################################

#10Kx10K

source("library_CIS_imagingData_08252020.R")
X_train <- read.csv("data_train_10X.csv", row.names=1)

dim(X_train)
#[1] 33538  6292
X_train <- t(X_train)
dim(X_train)
#[1]  6292 33538

y_train <- read.csv("y_Train_10X.csv")
dim(y_train)
y_train <- as.vector(y_train[,1])
length(y_train)
#[1] 6292

X_test <- read.csv("data_test_10X.csv", row.names=1)
dim(X_test)
#[1] 33538  1573
X_test <- t(X_test)
dim(X_test)
#[1]  1573 33538

table(y_train)

y_train_tmp <- rep(1, length(y_train))
y_train_tmp[which(y_train=="CD14+ monocytes")] <- 2
y_train_tmp[which(y_train=="CD16+ monocytes")] <- 3
y_train_tmp[which(y_train=="CD4+ T cells")] <- 4
y_train_tmp[which(y_train=="CD8+ T cells")] <- 5
y_train_tmp[which(y_train=="Dendritic cells")] <- 6
y_train_tmp[which(y_train=="NK cells")] <- 7
y_train_tmp[which(y_train=="unknown")] <- 8

table(y_train_tmp)

X_train[which(X_train==0, arr.ind=TRUE)] <- 0.01
logX <- log(X_train)

y0 <- as.vector(y_train_tmp)
X0 <- as.matrix(logX)
msd <- apply(X0, 2, sd) #### sd?
length(which(msd==0))
#[1] 12999
X0 <- X0[,-which(msd==0)]
X.new <- X0

dim(X0)
dim(X.new)
length(y0)

#try1  <- mLDA.Kclass(X0, y0, X.new, K=8, tau=20, alpha=0.7)
################## As the limited power of personal desktop, randomly select 2000 genes 
################## out of the dataset 
#X0.1 <- X0[,c(1:5000)]
#X0.2 <- X0[,c(5001:10000)]
#X0.3 <- X0[,c(10001:15000)]
#X0.4 <- X0[,c(15001:dim(X0)[2])]

#X.new.1 <- X.new[,c(1:5000)]
#X.new.2 <- X.new[,c(5001:10000)]
#X.new.3 <- X.new[,c(10001:15000)]
#X.new.4 <- X.new[,c(15001:dim(X0)[2])]

rs<-sample.int(ncol(X0), 2000)
#y0_rs <- y0
X0_rs <- X0[,as.factor(rs)]
X.new_rs <- X.new[,as.factor(rs)]

dim(X0_rs)
dim(X.new_rs)
length(y0_rs)

X0.1 <- X0_rs[,c(1:500)]
X0.2 <- X0_rs[,c(501:1000)]
X0.3 <- X0_rs[,c(1001:1500)]
X0.4 <- X0_rs[,c(1501:dim(X0_rs)[2])]

X.new.1 <- X.new_rs[,c(1:500)]
X.new.2 <- X.new_rs[,c(501:1000)]
X.new.3 <- X.new_rs[,c(1001:1500)]
X.new.4 <- X.new_rs[,c(1501:dim(X0_rs)[2])]

try1  <- mLDA.Kclass(X0_rs, y0, X.new_rs, K=8, tau=20, alpha=0.7)
length(which(try1$PredClass != y0_rs))
#[1] 1970

try1.1  <- mLDA.Kclass(X0.1, y0, X.new.1, K=8, tau=20, alpha=0.5)
length(which(try1.1$PredClass != y0))
#[1] 2413

All_set_1 <- try1.1$screenset
#screenset?
All_vox_1 <- colnames(X0.1[,All_set_1])
write.table(All_vox_1, file=paste("All_vox_1.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

length(All_vox_1)
#[1] 65

system.time(try1.2  <- mLDA.Kclass(X0.2, y0, X.new.2, K=8, tau=20, alpha=0.5))
##
length(which(try1.2$PredClass != y0))
#[1] 2354

All_set_2 <- try1.2$screenset
All_vox_2 <- colnames(X0.2[,All_set_2])
write.table(All_vox_2, file=paste("All_vox_2.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

length(All_vox_2)
#[1] 68

system.time(try1.3  <- mLDA.Kclass(X0.3, y0, X.new.3, K=8, tau=20, alpha=0.5))

length(which(try1.3$PredClass != y0))
#[1] 2535

All_set_3 <- try1.3$screenset
All_vox_3 <- colnames(X0.3[,All_set_3])
write.table(All_vox_3, file=paste("All_vox_3.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

length(All_vox_3)
#[1] 62

system.time(try1.4  <- mLDA.Kclass(X0.4, y0, X.new.4, K=8, tau=20, alpha=0.5))

length(which(try1.4$PredClass != y0))
#[1] 2457

All_set_4 <- try1.4$screenset
All_vox_4 <- colnames(X0.4[,All_set_4])
write.table(All_vox_4, file=paste("All_vox_4.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

length(All_vox_4)
#[1] 63

####################


All_sele_pred <-  rep(0, dim(X0_rs)[2])

for(r in 1:4){
  
  All_vox <- read.table(file=paste("All_vox_", r, ".txt", sep=""), header=FALSE)
  All_vox <- as.vector(All_vox[,1])
  all_id <- match(All_vox, colnames(X0_rs))
  All_sele_pred[all_id] <- All_sele_pred[all_id]+1
}


####################
### Prediction on training set

Fisher.M <- matrix(0, dim(X0_rs)[1], 8)

all_id <- which(All_sele_pred!=0)

X_All_train <- X0_rs[,all_id]
X_All_new <- X.new_rs[,all_id]

mean1.all <- apply(X_All_train[which(y0==1),], 2, mean) 
mean2.all <- apply(X_All_train[which(y0==2),], 2, mean) 
mean3.all <- apply(X_All_train[which(y0==3),], 2, mean)
mean4.all <- apply(X_All_train[which(y0==4),], 2, mean) 
mean5.all <- apply(X_All_train[which(y0==5),], 2, mean) 
mean6.all <- apply(X_All_train[which(y0==6),], 2, mean) 
mean7.all <- apply(X_All_train[which(y0==7),], 2, mean) 
mean8.all <- apply(X_All_train[which(y0==8),], 2, mean) 


mean1.all.m <- matrix(rep(mean1.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean1.all), byrow=TRUE)
mean2.all.m <- matrix(rep(mean2.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean2.all), byrow=TRUE)
mean3.all.m <- matrix(rep(mean3.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean3.all), byrow=TRUE)
mean4.all.m <- matrix(rep(mean4.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean4.all), byrow=TRUE)
mean5.all.m <- matrix(rep(mean5.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean5.all), byrow=TRUE)
mean6.all.m <- matrix(rep(mean6.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean6.all), byrow=TRUE)
mean7.all.m <- matrix(rep(mean7.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean7.all), byrow=TRUE)
mean8.all.m <- matrix(rep(mean8.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean8.all), byrow=TRUE)

Omega.r <- my.inv(cov(X_All_train))

Fisher.M[,1] <- as.matrix(X_All_new-mean1.all.m/2)%*%Omega.r%*%mean1.all
Fisher.M[,2] <- as.matrix(X_All_new-mean2.all.m/2)%*%Omega.r%*%mean2.all
Fisher.M[,3] <- as.matrix(X_All_new-mean3.all.m/2)%*%Omega.r%*%mean3.all
Fisher.M[,4] <- as.matrix(X_All_new-mean4.all.m/2)%*%Omega.r%*%mean4.all
Fisher.M[,5] <- as.matrix(X_All_new-mean5.all.m/2)%*%Omega.r%*%mean5.all
Fisher.M[,6] <- as.matrix(X_All_new-mean6.all.m/2)%*%Omega.r%*%mean6.all
Fisher.M[,7] <- as.matrix(X_All_new-mean7.all.m/2)%*%Omega.r%*%mean7.all
Fisher.M[,8] <- as.matrix(X_All_new-mean8.all.m/2)%*%Omega.r%*%mean8.all




predClass <- apply(Fisher.M, 1, which.max)
length(which(predClass!=y0))
#[1] 1649

### overall error rate 879/6292=0.14

length(which(predClass!=1 & y0==1))
#[1] 479
length(which(predClass!=2 & y0==2))
#[1] 1265
length(which(predClass!=3 & y0==3))
#[1] 4
length(which(predClass!=4 & y0==4))
#[1] 1932
length(which(predClass!=5 & y0==5))
#[1] 896
length(which(predClass!=6 & y0==6))
#[1] 26
length(which(predClass!=7 & y0==7))
#[1] 642
length(which(predClass!=8 & y0==8))
#[1] 801


####################
### Prediction on test set


X_test[which(X_test==0, arr.ind=TRUE)] <- 0.01
logXtest <- log(X_test)

X0test <- as.matrix(logXtest)
X0test <- X0test[,-which(msd==0)]
X.new <- X0test

dim(X.new)

X_All_new <- X.new_rs[,all_id]

Fisher.M <- matrix(0, dim(X.new_rs)[1], 8)

mean1.all.m <- matrix(rep(mean1.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean1.all), byrow=TRUE)
mean2.all.m <- matrix(rep(mean2.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean2.all), byrow=TRUE)
mean3.all.m <- matrix(rep(mean3.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean3.all), byrow=TRUE)
mean4.all.m <- matrix(rep(mean4.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean4.all), byrow=TRUE)
mean5.all.m <- matrix(rep(mean5.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean5.all), byrow=TRUE)
mean6.all.m <- matrix(rep(mean6.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean6.all), byrow=TRUE)
mean7.all.m <- matrix(rep(mean7.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean7.all), byrow=TRUE)
mean8.all.m <- matrix(rep(mean8.all, dim(X_All_new)[1]), dim(X_All_new)[1], length(mean8.all), byrow=TRUE)

Fisher.M[,1] <- as.matrix(X_All_new-mean1.all.m/2)%*%Omega.r%*%mean1.all
Fisher.M[,2] <- as.matrix(X_All_new-mean2.all.m/2)%*%Omega.r%*%mean2.all
Fisher.M[,3] <- as.matrix(X_All_new-mean3.all.m/2)%*%Omega.r%*%mean3.all
Fisher.M[,4] <- as.matrix(X_All_new-mean4.all.m/2)%*%Omega.r%*%mean4.all
Fisher.M[,5] <- as.matrix(X_All_new-mean5.all.m/2)%*%Omega.r%*%mean5.all
Fisher.M[,6] <- as.matrix(X_All_new-mean6.all.m/2)%*%Omega.r%*%mean6.all
Fisher.M[,7] <- as.matrix(X_All_new-mean7.all.m/2)%*%Omega.r%*%mean7.all
Fisher.M[,8] <- as.matrix(X_All_new-mean8.all.m/2)%*%Omega.r%*%mean8.all

predClass <- apply(Fisher.M, 1, which.max)

write.table(predClass, file="10Kx10K_Predict_results.txt", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
