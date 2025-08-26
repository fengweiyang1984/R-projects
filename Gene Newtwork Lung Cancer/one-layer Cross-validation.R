##### one-layer Cross-validation / tau=10 mu=100 alpha=0.5#####

#corss-validation
n=ncol(methylation_data1)
tau_r=10
nu_c=100
sepn = 5
lowb_upb = c(1:246)
lowb_upb1 = sample(lowb_upb, 50, replace=FALSE)
lowb_upb <- lowb_upb[!lowb_upb %in% lowb_upb1]

lowb_upb2 = sample(lowb_upb, 50, replace=FALSE)
lowb_upb <- lowb_upb[!lowb_upb %in% lowb_upb2]

lowb_upb3 = sample(lowb_upb, 50, replace=FALSE)
lowb_upb <- lowb_upb[!lowb_upb %in% lowb_upb3]

lowb_upb4 = sample(lowb_upb, 50, replace=FALSE)
lowb_upb <- lowb_upb[!lowb_upb %in% lowb_upb4]

lowb_upb5 = lowb_upb

lowb_upb = lowb_upb1

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
  
  data_set1_1 = data_set1[, lowb_upb]
  data_set2_1 = data_set2[, lowb_upb]
  dat_Test = cbind(data_set1_1,data_set2_1) 
  data_set1_1 = data_set1[, -lowb_upb]
  data_set2_1 = data_set2[, -lowb_upb]
  dat_Train = cbind(data_set1_1,data_set2_1)
  y_Train = as.vector(t(dat_Train[1,]))
  X_train = t(dat_Train[-1,])
  y_Test = as.vector(t(dat_Test[1,]))
  X_test = t(dat_Test[-1,])
  
  #seperate Z
  data_set1 = Z[,colnames(Z)%in%name1$sample]
  data_set2 = Z[,colnames(Z)%in%name2$sample]
  data_set1_1 = data_set1[, lowb_upb]
  data_set2_1 = data_set2[, lowb_upb]
  dat_Test = cbind(data_set1_1,data_set2_1) 
  data_set1_1 = data_set1[, -lowb_upb]
  data_set2_1 = data_set2[, -lowb_upb]
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

predClass <- apply(Fisher.M, 1, which.max)
test_error=length(which(predClass!=y))
test_error

#test_error_1=1  
#test_error_2=1 
#test_error_2=1 
#test_error_2=1
#test_error_2=2
(1+1+1+1+2)/492

(1+0+1+0+1)/492

(1+1+1+0+1)/492
