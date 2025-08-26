##################################################################################################
##### case-control classification analysis
##################################################################################################
##################################################################################################
setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Real data example")

##########Extract the known genes in disease related KEGG pathways
library("KEGGREST")
namestrue <- as.data.frame(NULL)
#Get the list of numbers, gene symbols and gene description
names <- keggGet("hsa04151")[[1]]$GENE
#Delete the gene number by deleting every other line
namesodd <-  names[seq(0,length(names),2)]
#Create a substring deleting everything after the ; on each line (this deletes the gene description).
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

names <- keggGet("hsa04919")[[1]]$GENE
namesodd <-  names[seq(0,length(names),2)]
temp <- as.data.frame(gsub("\\;.*","",namesodd))
namestrue <- rbind(namestrue,temp)

#names <- keggGet("hsa05206")[[1]]$GENE
#namesodd <-  names[seq(0,length(names),2)]
#temp <- as.data.frame(gsub("\\;.*","",namesodd))
#namestrue <- rbind(namestrue,temp)

#names <- keggGet("hsa05200")[[1]]$GENE
#namesodd <-  names[seq(0,length(names),2)]
#temp <- as.data.frame(gsub("\\;.*","",namesodd))
#namestrue <- rbind(namestrue,temp)

temp <- unique(toupper(namestrue[,1]))
#subset X matrix based on gene symbols from namestrue
X_input <- exprALL
X_input <- X_input[,colnames(X_input)%in%temp]
X_input <- apply(X_input, 2, as.numeric)
X_input <- as.data.frame(X_input);
#X_input <- apply(X_input, 2, scale)
X_input <- round(X_input,3)
X_input <-apply(X_input,2,function(x)((x-mean(x))/sd(x)))
X_input <- X_input[ , colSums(is.na(X_input))==0]
#X_input <- X_input[,-(which(colSums(X_input)==0))]
X_input <- X_input[,-which(colnames(X_input) %in% c("IFNA4"))]
rownames(X_input) <- rownames(exprALL)

#dim(X_input); str(X_input)
y_input <- as.matrix(cbind(survival_OV[,4],survival_OV[,3])); 
colnames(y_input) <- c("time","status"); rownames(y_input)<-survival_OV[,1]
idx <- match(rownames(X_input),rownames(y_input))
y_input <- y_input[idx,]; 
y_input[which(y_input[,1] ==0 )]<-0.1
y_input[is.na(y_input)] <- 0

#####cross validation to determine the best lambda in group lasso
IRLS.IP.gglasso.crossvalidation <-
  function(X2, y2, n, sepn =10,lambda_max = 0.1){
    valid_result <- NULL
    for (k in seq(0.06, lambda_max, length.out = 7)) {
    valid_output = NULL
    for (tt in 1:sepn) {
      #### Preparing the training & testing dataset
      lowb1 = ceiling((n)/sepn*(tt-1))+1
      upb1 = min(ceiling((n)/sepn*tt),n)

      X_train = X2[-lowb1:-upb1,]
      X_test = X2[lowb1:upb1,] 
      
      y_Train = y2[-lowb1:-upb1,]
      y_Test =  y2[lowb1:upb1,]
      
      XtY <- cor(X_train,y_Train[,1])
      XtY[which(abs(XtY)<0.01)]<-0
      XtY <- as.data.frame(XtY)
      XtX <- cor(X_train)
      Sigma <- cov(X_train)
      Sigma <- round(Sigma,digits = 3)
      library(matlib)
      Sigma_inv <- solve(Sigma)
      #dim(Sigma_inv)
      Sigma_inv[which(abs(XtX)<=0.5)] <- 0
      event_raw1  <-as.data.frame(y_Train[,2]); colnames(event_raw1) <- "event"
      time_raw1 <- as.numeric(y_Train[,1])
      #X_train <-apply(X_train,2,function(x)((x-mean(x))/sd(x)))
      rownames(X_train) <- rownames(X2[-lowb1:-upb1,])
      
      #fit model
      cox_mine <- IRLS_IP_gglasso(x=X_train,times=time_raw1,event=event_raw1,lambda_lasso=k)
      #predict and calculate C-index
      event_raw1  <-as.data.frame(y_Test[,2]); event_raw1 <- sapply(event_raw1,as.numeric);colnames(event_raw1) <- "event"
      time_raw1 <- as.numeric(y_Test[,1])
      cindex_temp <- cindex_calcu(time2=time_raw1, status2=event_raw1, X2=X_test, B2=cox_mine)[[3]]
      
      valid_output <- rbind(valid_output,cindex_temp)
    }
  temp1 <- as.vector(c(k,mean(valid_output)))
  valid_result <- rbind(valid_result,temp1)
  print(k)
  }
  return(valid_result)
}
source("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Code/Liabary/Proportional hazards model V3_1.R")
valid_result <- IRLS.IP.gglasso.crossvalidation(X2=X_input, y2=y_input, n=nrow(X_input))
write.csv(valid_result,"valid_result.csv")

library(glmnet)
temp.lasso.cvfit <- cv.glmnet(x=X_input,y=y_input,family = "cox", type.measure = "C",alpha = 1,nfolds = 10)
plot(temp.lasso.cvfit)
temp.lasso.cvfit$lambda.1se ###or temp.lasso.cvfit$lambda.min ### or
#> temp.lasso.cvfit$lambda.1se
#[1] 0.02867632
temp.EN.cvfit <- cv.glmnet(x=X_input,y=y_input, family = "cox", type.measure = "C",alpha = 0.5,nfolds = 10)
plot(temp.EN.cvfit)
temp.EN.cvfit$lambda.1se ###or temp.EN.cvfit$lambda.min ### or 
#> temp.lasso.cvfit$lambda.1se
#[1] 0.06008352
temp.Ri.cvfit <- cv.glmnet(x=X_input,y=y_input, family = "cox", type.measure = "C",alpha = 0,nfolds = 10)
plot(temp.Ri.cvfit)
temp.Ri.cvfit$lambda.1se ###or temp.Ri.cvfit$lambda.min ### or 
#[1] 7.795909


#correlation betweeen X and y
XtY <- cor(X_input,y_input)
XtY[which(abs(XtY)<=0.5)]<-0
XtY <- as.data.frame(XtY)

#correlation of X input
XtX <- pre_network
#covariance of X
Sigma <- cov(X_input)
Sigma <- round(Sigma,digits = 5)
library(matlib)
#precision matrix
Sigma_inv <- solve(Sigma)
dim(Sigma_inv)
#use a threshold on correlation of X to threshold precision matrix
Sigma_inv[which(XtX==0)] <- 0
event_raw1  <-as.data.frame(y_input[,2]); colnames(event_raw1) <- "event"
time_raw1 <- as.numeric(y_input[,1])

#model fitting
#HDnetCox
source("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 2/Code/Liabary/Proportional hazards model V3_1.R")
try(Cox_mine <- IRLS_IP_gglasso(x=X_input,times=time_raw1,event=event_raw1,lambda_lasso=0.04))
length(which(Cox_mine!=0))

library("glmnet") 
#lasso
lamb <- temp.lasso.cvfit$lambda.1se ### or lamb <- temp.lasso.cvfit$lambda.min
temp.lasso <- glmnet(x=X_input, y=y_input, family = "cox",lambda=lamb,alpha = 1)
cox_laaso <- as.matrix(temp.lasso[["beta"]])
#length(which(cox_laaso!=0))
#[1] 21

#elastic net
lamb <- temp.EN.cvfit$lambda.1se ### or lamb <- temp.lasso.cvfit$lambda.min
temp.en <- glmnet(x=X_input, y=y_input, family = "cox",lambda=lamb,alpha = 0.5)
cox_EN <- as.matrix(temp.en[["beta"]])
#length(which(cox_EN!=0))

###
temp <- cbind(glm_mine,glm_mine1,cox_laaso,cox_EN) 
colnames(temp) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_EN")
beta_estimation <- temp
y_prediction <- NULL
for (i in 1:ncol(beta_estimation)){
  temp <- exp(X_input%*%beta_estimation[,i])/(1+exp(X_input%*%beta_estimation[,i]))
  y_prediction <- cbind(y_prediction,temp)
}
colnames(y_prediction) <- c("IRLS_IP_gglasso","IRLS_IP_sparsegl","glm_laaso","glm_EN")
y_prediction <- as.data.frame(y_prediction)
y_prediction$true_outcome <- caseControlResponse





###################   Visualization of pathways
library(pathview)
data(gse16873.d)
data(demo.paths)
i <- 1
#rownames are hsa id for genes, first column is prerank score.
pathway_comparison <- as.matrix(temp)


v.out <- pathview(gene.data = pathway_comparison[,1], pathway.id = "04151",species = "hsa", 
                  out.suffix = "Cox_PI3K-Akt", kegg.native = T,same.layer = F,key.pos = "bottomright",  sign.pos = "topright",
                  low = list(gene = "gray80", cpd = "blue"), 
                  mid = list(gene = "gray58", cpd = "gray"), 
                  high = list(gene = "gray30", cpd = "yellow"))

#####survival kaplan meier based on the selected gene from HDnetCox ####################
library(survival)
library(ggplot2)
library(survminer)

y_input_KM <- as.data.frame(y_input)
X_input_KM <- as.data.frame(X_input)
#based on the gene expression to separate data to low/high classes
PredClass <- rep(0,dim(y_input_KM)[1])
#MAPK3, ERBB2, CDKN1B, and AKT1
PredClass[which(X_input_KM$AKT1 >mean(X_input_KM$AKT1))] <- 1
y_input_KM$OSD <- PredClass
#KM plot
km_fit <- survfit(Surv(time,status) ~ OSD, data=y_input_KM)
ggsurv <-ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Overall survival probability", risk.table = TRUE, tables.theme = clean_theme(),
  linetype = c("solid","dashed"), palette = c("red", "lightgreen"))

ggsurv$table <- ggsurv$table +
  scale_y_discrete(breaks = c("OSD=0", "OSD=1"), 
                   labels = rev(c("Low", "High") ))
ggsurv

km_fit <- survfit(Surv(time,status) ~ 1, data=y_input_KM)
ggsurv <-ggsurvplot(
  fit = km_fit, 
  xlab = "Survival Days", 
  ylab = "Overall survival probability", risk.table = TRUE, tables.theme = clean_theme(),
  linetype = "solid", palette = c("#02080c"))

ggsurv
