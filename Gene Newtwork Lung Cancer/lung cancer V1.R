getwd()
setwd("C:\\Users\\Fengwei Yang\\Downloads")
memory.limit(size = 160000)
memory.size(max=T)
require(ggplot2)

#input survival data
TCGA_LUAD_survival <- read.delim("TCGA-LUAD.survival.tsv.gz", quote="")
write.table(TCGA_LUAD_survival, "survival.txt", append = FALSE, sep = " ", dec = ".",
            row.names = TRUE, col.names = TRUE)

TCGA.LUAD.GDC_phenotype <- read.delim("TCGA-LUAD.GDC_phenotype.tsv.gz")
TCGA_LUAD_methylation450 <- read.table(file="TCGA-LUAD.methylation450.tsv.gz",sep = "\t", header = TRUE)
rown = TCGA_LUAD_methylation450[,1]
#abc<-TCGA_LUAD_methylation450[1,]
methylation = TCGA_LUAD_methylation450[,-1]
methylation = as.data.frame(methylation)
methylation = as.matrix(methylation)

#extract the patient in survival and methlation data
patientid_m <- as.data.frame(colnames(methylation))

#sur_data <- t1
summary(TCGA_LUAD_survival$OS.time)
Z_1 = as.data.frame(t(Z))
summary(Z_1$number_pack_years_smoked)
summary(Z_1$age_at_initial_pathologic_diagnosis)
boxplot(TCGA_LUAD_survival$OS.time)
ggplot(TCGA_LUAD_survival, aes(y=OS.time))+
  geom_boxplot(outlier.colour=NA) #+ coord_cartesian(ylim = c(0, 2000))
########################
######################## The mean of OS.time is 900.9, and median is 657.5
########################
methylation_data <- TCGA_LUAD_methylation450[,c(-29,-91,-106,-132,-153,-163,-216,-288,-298,-373,-440)]
rownames(methylation_data) = rown[,1]

#check the distribution of methylation data
#remove the row with NA 
for (i in 1:ncol(methylation_data)) {
  methylation_data1 <-subset(methylation_data, !is.na(methylation_data[,i]))  
}
#pick up first 10k rows from dataset
illuminaMethyl450_hg38_GDC <- read.table("illuminaMethyl450_hg38_GDC")
illuminaMethyl450_hg38_GDC = illuminaMethyl450_hg38_GDC[,c(1,3)]
colnames(illuminaMethyl450_hg38_GDC) =  c("Marker","Chromsome")
cc=22
chr_n=paste("chr", cc, sep="")
marker_chr = illuminaMethyl450_hg38_GDC[illuminaMethyl450_hg38_GDC$Chromsome==chr_n,]
methylation_data_10k = methylation_data1[rownames(methylation_data1)%in%marker_chr$Marker,]

# log tranformation 
methylation_data_10k[which(methylation_data_10k==0, arr.ind=TRUE)] <- 0.01
         # Correct the cell whose value is NA/na/N/a
methylation_data_10k[is.na(methylation_data_10k)]=0.01
methylation_data_10k[which(methylation_data_10k=="na", arr.ind=TRUE)]=0.01
methylation_data_10k[which(methylation_data_10k=="N/A", arr.ind=TRUE)]=0.01
log_methylation_data_10k <- log(methylation_data_10k)
log_methylation_data_10k = t(scale(t(log_methylation_data_10k), center = TRUE, scale = TRUE))
#abc = log_methylation_data_10k[1:2,]
#rownames(abc)= c("samplen","OSD")
#abc[1,]=colnames(abc)
#abc= as.data.frame(t(abc))
#write.csv(abc, file="abc.csv")
#write.csv(TCGA.LUAD.survival, file="TCGA.LUAD.survival.csv")
#abc = as.data.frame(t(abc))
#colnames(abc)=c("samplen","OSD")
#colnames(abc) = as.factor(abc[1,])
log_methylation_data_10k = rbind(abc,log_methylation_data_10k)

#seperate training & testing dataset
name1 = TCGA.LUAD.survival[TCGA.LUAD.survival$OSD==1,]
name2 = TCGA.LUAD.survival[TCGA.LUAD.survival$OSD==0,]
data_set1 = log_methylation_data_10k[,colnames(log_methylation_data_10k)%in%name1$sample]
data_set2 = log_methylation_data_10k[,colnames(log_methylation_data_10k)%in%name2$sample]
data_set2[1,] = 2
#data_set1[is.na(data_set1)]=log(0.01)
#data_set2[is.na(data_set2)]=log(0.01)
#rm(log_methylation_data_10k)

##### double-layers Cross-validation #####
#require(dplyr)
n=ncol(methylation_data_10k)
tau_seq = 5
nu_seq = 5
sepn = 6
mLDA.pair.crossvalidation <-
  function(data_set1, data_set2, n, sepn = sepn,tau_seq = tau_seq,nu_seq = nu_seq){
    output = data.frame()
    temp_output = data.frame()
    temp_output = data.frame(vi=0,tau=0, nu=0, error =0)
    temp_output1 = data.frame()
    valid_output = data.frame()
    for (i in 1:sepn) {
      #### Preparing the training & testing dataset
      #i=2
      lowb = n/(sepn*2)*(i-1)+1
      upb = n/(sepn*2)*i
      data_set1_1 = data_set1[, lowb:upb]
      data_set2_1 = data_set2[, lowb:upb]
      dat_Test = cbind(data_set1_1,data_set2_1) 
      data_set1_1 = data_set1[, -lowb:-upb]
      data_set2_1 = data_set2[, -lowb:-upb]
      dat_Train = cbind(data_set1_1,data_set2_1)
      y_Train = as.vector(t(dat_Train[1,]))
      X_train = t(dat_Train[-1,])
      y_Test = as.vector(t(dat_Test[1,]))
      X_test = t(dat_Test[-1,])
      #### mLDA fit ####
      source("library_CIS_imagingData_0820.R")
      for (r in 1:tau_seq) {
        #r=2
        tau_r = 5 * r 
        for (c in 1:nu_seq) {
          #c=4
          nu_c = 50 * c
          try1  <- mLDA.pair(X_train, y_Train, X_test, tau=tau_r, nu=nu_c)
          temp_output[1] = i
          temp_output[2] = tau_r  
          temp_output[3] = nu_c  
          temp_output[4] = length(which(try1$PredClass != y_Test))  
          temp_output1 = rbind(temp_output1,temp_output)
        }
      }
      #temp_output1 = subset(temp_output1,temp_output1$error==with(temp_output1, min(temp_output1$error)))
    }
    valid_output = temp_output1
    #temp_output1 = subset(valid_output,valid_output$error==with(valid_output, min(valid_output$error)))
    #i = temp_output1[1]
    #lowb = n/(sepn*2)*(i-1)+1
    #upb = n/(sepn*2)*i
    #data_set1_1 = data_set1[, lowb:upb]
    #data_set2_1 = data_set2[, lowb:upb]
    #dat_Test = cbind(data_set1_1,data_set2_1) 
    #data_set1_1 = data_set1[, -lowb:-upb]
    #data_set2_1 = data_set2[, -lowb:-upb]
    #dat_Train = cbind(data_set1_1,data_set2_1)
    #y_Train = as.vector(t(dat_Train[1,]))
    #X_train = t(dat_Train[-1,])
    #y_Test = as.vector(t(dat_Test[1,]))
    #X_test = t(dat_Test[-1,])
    #tau_r = temp_output[2]   
    #nu_c = temp_output[3] 
    #try1  <- mLDA.pair(X_train, y_Train, X_test, tau=tau_r, nu=nu_c)
    return(valid_output)
    valid_output
  }
mLDA.pair.crossvalidation(data_set1=data_set1, data_set2=data_set2, n=ncol(methylation_data_10k), sepn = sepn,tau_seq = tau_seq,nu_seq = nu_seq)
min(valid_output$error)
subset(valid_output,valid_output$error==min(valid_output$error))
write.csv(valid_output,file=paste("chr", cc,".csv", sep=""))

##### one-layer Cross-validation / tau=10 mu=100 #####

#corss-validation
n=ncol(methylation_data_10k)
tau_r=10
nu_c=100
sepn = 5
i=1
lowb1 = floor(n/(sepn*2)*(i-1)+1)+10
upb1 = ceiling(n/(sepn*2)*i)+10
i=2
lowb2 = upb1+1
upb2 = ceiling(n/(sepn*2)*i)+11
i=3
lowb3 = upb2+1
upb3 = ceiling(n/(sepn*2)*i)+12
i=4
lowb4 = upb3+1
upb4 = ceiling(n/(sepn*2)*i)+13
i=5
lowb5 = upb4+1
upb5 = ceiling(n/(sepn*2)*i)

lowb = lowb1
upb = upb1

for (cc in 1:22) {
  #cc=20
  chr_n=paste("chr", cc, sep="")
  marker_chr = illuminaMethyl450_hg38_GDC[illuminaMethyl450_hg38_GDC$Chromsome==chr_n,]
  methylation_data_10k = methylation_data1[rownames(methylation_data1)%in%marker_chr$Marker,]
  
  # Correct the cell whose value is NA/na/N/a
  methylation_data_10k[which(methylation_data_10k==0, arr.ind=TRUE)] <- 0.01
  methylation_data_10k[is.na(methylation_data_10k)]=0.01
  methylation_data_10k[which(methylation_data_10k=="na", arr.ind=TRUE)]=0.01
  methylation_data_10k[which(methylation_data_10k=="N/A", arr.ind=TRUE)]=0.01
  log_methylation_data_10k <- log(methylation_data_10k)
  # scale dataset 
  log_methylation_data_10k = t(scale(t(log_methylation_data_10k), center = TRUE, scale = TRUE))
  #log_methylation_data_10k = rbind(abc,log_methylation_data_10k)
  log_methylation_data_10k = rbind(Z[1,],log_methylation_data_10k)

  #seperate training & testing dataset
  name1 = TCGA.LUAD.survival[TCGA.LUAD.survival$OSD==1,]
  name2 = TCGA.LUAD.survival[TCGA.LUAD.survival$OSD==0,]
  data_set1 = log_methylation_data_10k[,colnames(log_methylation_data_10k)%in%name1$sample]
  data_set2 = log_methylation_data_10k[,colnames(log_methylation_data_10k)%in%name2$sample]
  data_set2[1,] = 2
  
  data_set1_1 = data_set1[, lowb:upb]
  data_set2_1 = data_set2[, lowb:upb]
  dat_Test = cbind(data_set1_1,data_set2_1) 
  data_set1_1 = data_set1[, -lowb:-upb]
  data_set2_1 = data_set2[, -lowb:-upb]
  dat_Train = cbind(data_set1_1,data_set2_1)
  y_Train = as.vector(t(dat_Train[1,]))
  X_train = t(dat_Train[-1,])
  y_Test = as.vector(t(dat_Test[1,]))
  X_test = t(dat_Test[-1,])
  
  #seperate Z
  data_set1 = Z[,colnames(Z)%in%name1$sample]
  data_set2 = Z[,colnames(Z)%in%name2$sample]
  data_set1_1 = data_set1[, lowb:upb]
  data_set2_1 = data_set2[, lowb:upb]
  dat_Test = cbind(data_set1_1,data_set2_1) 
  data_set1_1 = data_set1[, -lowb:-upb]
  data_set2_1 = data_set2[, -lowb:-upb]
  dat_Train = cbind(data_set1_1,data_set2_1)
  Z_train = t(dat_Train[-1,])
  Z_test = t(dat_Test[-1,])
  
  #### mLDA fit ####
  source("library_CIS_imagingData_0820.R")
  #try1  <- mLDA.pair(X_train, y_Train, X_test, Z=Z_train, Z.new=Z_test, tau=tau_r, nu=nu_c)
  try1  <-mLDA.pair(X=X_train, y=y_Train, X.new=X_test, Z=Z_train, Z.new=Z_test,  tau=tau_r, alpha=0.5, nu=nu_c)
  saveRDS(try1,file=paste("MI_CpG_chr", cc, ".rds", sep=""))
  }

               #Combine ALL_set and compute the test error
##### considering number_pack_years_smoked & age_at_initial_pathologic_diagnosis as covariances #####
#temp = data.frame()
#temp = TCGA.LUAD.GDC_phenotype$age_at_initial_pathologic_diagnosis
#Z = as.data.frame(TCGA.LUAD.GDC_phenotype[,1])
#Z$number_pack_years_smoked = TCGA.LUAD.GDC_phenotype$number_pack_years_smoked 
#Z$age_at_initial_pathologic_diagnosis = TCGA.LUAD.GDC_phenotype$age_at_initial_pathologic_diagnosis
#colnames(Z)[1]="sample"
#write.csv(Z,"Z.csv")
Z[is.na(Z)]=0
Z_1=Z[,1]
Z_2=Z[,-1]
Z_2 = scale(Z_2, center = TRUE, scale = TRUE)
Z = t(cbind(Z_1,Z_2))
rownames(Z)[1]="OSD"

#marker_all = NULL
Fisher_combined0 = NULL
Fisher_combined1 = NULL

for(cc in 1:22){
  ###############Fisher.M <- matrix(0, dim(X_test)[2], 2)
  #cc=2
  chr_n=paste("chr", cc, sep="")
  try1 = readRDS(file=paste("MI_CpG_chr", cc, ".rds", sep=""))
  marker_chr = illuminaMethyl450_hg38_GDC[illuminaMethyl450_hg38_GDC$Chromsome==chr_n,]
  methylation_data_10k = methylation_data1[rownames(methylation_data1)%in%marker_chr$Marker,]
  
  # Correct the cell whose value is NA/na/N/a
  methylation_data_10k[which(methylation_data_10k==0, arr.ind=TRUE)] <- 0.01
  methylation_data_10k[is.na(methylation_data_10k)]=0.01
  methylation_data_10k[which(methylation_data_10k=="na", arr.ind=TRUE)]=0.01
  log_methylation_data_10k <- log(methylation_data_10k)
  # scale dataset 
  log_methylation_data_10k = t(scale(t(log_methylation_data_10k), center = TRUE, scale = TRUE))
  All_CpG = rownames(as.data.frame(try1$network.nodes))
  X_All = log_methylation_data_10k[rownames(log_methylation_data_10k)%in%All_CpG,]
  X_All_1 = rbind(Z[1,],X_All)
  
  data_set1 = X_All_1[,colnames(X_All_1)%in%name1$sample]
  data_set2 = X_All_1[,colnames(X_All_1)%in%name2$sample]
  data_set2[1,] = 2
  
  data_set1_1 = data_set1[, -lowb:-upb]
  data_set2_1 = data_set2[, -lowb:-upb]
  dat_Train = cbind(data_set1_1,data_set2_1)
  data_set1_1 = data_set1[, lowb:upb]
  data_set2_1 = data_set2[, lowb:upb]
  dat_Test = cbind(data_set1_1,data_set2_1)
  X_test = dat_Test[-1,]
  Z_test = Z[,colnames(Z)%in%colnames(X_test)]
  y=dat_Test[1,]

  Fisher.M <- matrix(0, dim(X_test)[2], 2)
  temp0=dat_Test[,which(dat_Test[1,]==1)]
  temp0=temp0[-1,]
  mean0.all <- as.matrix(apply(temp0, 1, mean))
  temp1=dat_Test[,which(dat_Test[1,]==2)]
  temp1=temp1[-1,]
  mean1.all <- as.matrix(apply(temp1, 1, mean)) 
  
  mean0.all.m <- matrix(rep(mean0.all, dim(X_test)[2]), nrow=dim(X_test)[1], byrow=F)
  mean1.all.m <- matrix(rep(mean1.all, dim(X_test)[2]), nrow=dim(X_test)[1], byrow=F)
  
  Omega <- try1$network.Omega
  #Omega <- qr.solve(cov(X_test))
  Omega.mi <- matrix(0, dim(X_All)[1], dim(X_All)[1]) 
  var.mi <- apply(X_All, 1, var)
  diag(Omega.mi) <- var.mi
  #t_without_z = (t(X_test-mean0.all.m/2))%*%Omega%*%mean0.all
  
  Fisher.M[,1] <- Fisher.M[,1] + (t(X_test-mean0.all.m/2))%*%Omega%*%mean0.all
  Fisher.M[,2] <- Fisher.M[,2] + (t(X_test-mean1.all.m/2))%*%Omega%*%mean1.all

  #dim(mean0.all.m)
  Fisher_combined0 = as.data.frame(cbind(Fisher_combined0,Fisher.M[,1]))
  Fisher_combined1 = as.data.frame(cbind(Fisher_combined1,Fisher.M[,2]))
}
#write.table(Fisher_combined0, file="Fisher_combined0_f5.csv", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
#write.table(Fisher_combined1, file="Fisher_combined1_f5.csv", quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
Fisher.M[,1]=apply(Fisher_combined0,1,sum)
Fisher.M[,2]=apply(Fisher_combined1,1,sum)

temp0=Z_test[,which(Z_test[1,]==0)]
temp0=temp0[-1,]
mean0.Z <- as.matrix(apply(temp0, 1, mean))
temp1=Z_test[,which(Z_test[1,]==1)]
temp1=temp1[-1,]
mean1.Z <- as.matrix(apply(temp1, 1, mean))
Z_test =Z_test[-1,]

mean0.Z.m <- matrix(rep(mean0.Z, dim(Z_test)[2]), nrow=dim(Z_test)[1], byrow=F)
mean1.Z.m <- matrix(rep(mean1.Z, dim(Z_test)[2]), nrow=dim(Z_test)[1], byrow=F)

#ZOmega <- qr.solve(cov(t(Z_test)))
ZOmega <- cov(t(Z_test))

Fisher.M[,1] <- Fisher.M[,1] + (t(Z_test-mean0.Z.m/2))%*%ZOmega%*%mean0.Z
Fisher.M[,2] <- Fisher.M[,2] + (t(Z_test-mean1.Z.m/2))%*%ZOmega%*%mean1.Z

predClass <- apply(Fisher.M, 1, which.max)
test_error=length(which(predClass!=y))
test_error
#test_error_1=9  
#test_error_2=1 
#test_error_2=0 
#test_error_2=0
#test_error_2=0
(9+16+18+18+13)/492
#test_error_rate=0.1382114  / 0.1382114
y[1,]<-predClass
#predClass_all <- y
predClass_all <- cbind(predClass_all,y)
#Use Z to calculate in each chrom                
#test_error_1=9    
#test_error_2=0 
#test_error_2=18 
#test_error_2=18 
#test_error_2=13 
#test_error_rate=0.152439  

################################################################
################################################################
#Get the Cpg sites for whole dataset
for (cc in 1:22) {
  #cc=20
  chr_n=paste("chr", cc, sep="")
  marker_chr = illuminaMethyl450_hg38_GDC[illuminaMethyl450_hg38_GDC$Chromsome==chr_n,]
  methylation_data_10k = methylation_data1[rownames(methylation_data1)%in%marker_chr$Marker,]
  
  # Correct the cell whose value is NA/na/N/a
  methylation_data_10k[which(methylation_data_10k==0, arr.ind=TRUE)] <- 0.01
  methylation_data_10k[is.na(methylation_data_10k)]=0.01
  methylation_data_10k[which(methylation_data_10k=="na", arr.ind=TRUE)]=0.01
  methylation_data_10k[which(methylation_data_10k=="N/A", arr.ind=TRUE)]=0.01
  log_methylation_data_10k <- log(methylation_data_10k)
  # scale dataset 
  log_methylation_data_10k = t(scale(t(log_methylation_data_10k), center = TRUE, scale = TRUE))
  #log_methylation_data_10k = rbind(abc,log_methylation_data_10k)
  log_methylation_data_10k = rbind(abc,log_methylation_data_10k)
  
  name1 = TCGA.LUAD.survival[TCGA.LUAD.survival$OSD==1,]
  name2 = TCGA.LUAD.survival[TCGA.LUAD.survival$OSD==0,]
  data_set1 = log_methylation_data_10k[,colnames(log_methylation_data_10k)%in%name1$sample]
  data_set2 = log_methylation_data_10k[,colnames(log_methylation_data_10k)%in%name2$sample]
  data_set2[1,] = 2
  log_methylation_data_10k = cbind(data_set1,data_set2)
  
  #seperate training & testing dataset
  X_train = t(log_methylation_data_10k[-1,])
  y_Train = as.vector(t(log_methylation_data_10k[1,]))
  X_test = t(log_methylation_data_10k[-1,])
  Z_train = t(Z[-1,])
  Z_test = t(Z[-1,])
  
  #### mLDA fit ####
  source("library_CIS_imagingData_0820.R")
  #try1  <- mLDA.pair(X_train, y_Train, Z=Z_train,X_test, Z.new=Z_test, tau=tau_r, nu=nu_c)
  try1  <-mLDA.pair(X=X_train, y=y_Train, X.new=X_test, Z=Z_train, Z.new=Z_test,  tau=tau_r, alpha=0.5, nu=nu_c) 
  saveRDS(try1,file=paste("MI_CpG_chr", cc, ".rds", sep=""))
}
