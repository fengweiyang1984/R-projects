setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 3/Paper work/Figures")
#X_input <- Logistic_Data
#Logistic_Data <- regionM[1:184,]

#regionM <- NULL
#for (t in 1:length(Tb_mask_non0)){
#  tmpM <- All_image[, which(mask_non0==names(Tb_mask_non0)[t])]
#  tmp.v <- apply(tmpM, 1, mean)
#  regionM <- cbind(regionM, tmp.v)
#}

#All_image <- rbind(AD_image, MCI_image, NORM_image) 

#AD_image <- AD_b[which(!is.na(AD_mid)), ]
#MCI_image <- MCI_b[which(!is.na(MCI_mid)), ]

#AD_b <- (dat$AD)$bl
#MCI_b <- (dat$MCI)$bl
par(mfrow=c(2,2))

#X is right(+), left(-)
ImgX <- AD1.arr[46,,]
p1 <- MyHeatMapADNI(ImgX)

#Y is front(+) back(-)
ImgY <- AD1.arr[,54,]
p2 <- MyHeatMapADNI(ImgY)

#Z is top(+) down(-)
ImgZ <- AD1.arr[,,45]
p3 <- MyHeatMapADNI(ImgZ)

library(ggplot2)
library(gtable)
library(gridExtra)

#grid.arrange(p1,p2,p3,ncol=2)
#grid.arrange(P1, P2, P3, ncol = 2, nrow = 2, layout_matrix= rbind(c(1,2), 3))

#install.packages("ggseg3d")
library(ggseg3d)
library(dplyr)
library(tidyr)
vignette("ggseg3d")
data(aseg_3d)

##########################Glass brain figures#################################
##########################Glass brain figures#################################
##########################Glass brain figures#################################

######Brodmann’s area 24/32/33
######Brodmann area 24 is part of the anterior cingulate.
######Brodmann area 32, also known in the human brain as the dorsal anterior cingulate area 32
######Brodmann area 33, pregenual area 33, is a subdivision of the cytoarchitecturally defined cingulate region.
someData = dk_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>% 
  filter(label=="rh_caudalanteriorcingulate"|label=="rh_posteriorcingulate"|label=="rh_rostralanteriorcingulate"|
           label=="lh_caudalanteriorcingulate"|label=="lh_posteriorcingulate"|label=="lh_rostralanteriorcingulate") %>% 
  na.omit() %>% 
  mutate(p = seq(1, nrow(.)))

ggseg3d(.data = someData, 
        atlas = dk_3d,
        colour = "p", text = "p", 
        na.alpha= .2) %>% 
  add_glassbrain() %>% 
  remove_axes() 

######Brodmann’s area 23/31
######Brodmann area 23 (BA23) is a region in the brain that lies inside the posterior cingulate cortex. 
######Brodmann area 31, also known as dorsal posterior cingulate area 31, is a subdivision of the cytoarchitecturally defined cingulate region of the cerebral cortex
someData = dk_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>% 
  filter(label=="rh_isthmuscingulate"|
           label=="lh_isthmuscingulate") %>% 
  na.omit() %>% 
  mutate(p = seq(1, nrow(.)))

ggseg3d(.data = someData, 
        atlas = dk_3d,
        colour = "p", text = "p", 
        na.alpha= .2) %>% 
  add_glassbrain() %>% 
  remove_axes() 

######Brodmann’s area 35
######Brodmann area 35 (BA23) is known as perirhinal area 35. It is a subdivision of the cytoarchitecturally defined hippocampal region of the cerebral cortex. . 
somData_aseg = aseg_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>% 
  filter(label=="Left-Hippocampus"|label=="Right-Hippocampus") %>% 
  mutate(p = seq(1, nrow(.))) 

ggseg3d(.data = somData_aseg, atlas = aseg_3d,colour = "p", text = "p", 
        na.alpha= .2) %>% 
  add_glassbrain() %>% 
  remove_axes() %>% 
  pan_camera("right lateral")

######Brodmann’s area 39
######Brodmann area 39 is part of the parietal cortex in the human brain 
someData = dk_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>% 
  filter(label=="lh_inferiorparietal") %>% 
  na.omit() %>% 
  mutate(p = seq(1, nrow(.)))


ggseg3d(.data = someData, 
        atlas = dk_3d,
        colour = "p", text = "p",na.alpha= .2) %>% 
  add_glassbrain() %>% 
  remove_axes() 

######Putamen
somData_aseg = aseg_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>% 
  filter(label=="Left-Putamen"|label=="Right-Putamen") %>% 
  na.omit() %>% 
  mutate(p = seq(1, nrow(.))) 

ggseg3d(.data = somData_aseg, atlas = aseg_3d,colour = "p", text = "p", 
        na.alpha= .2) %>% 
  add_glassbrain() %>% 
  remove_axes() %>% 
  pan_camera("right lateral")


##########################Averaged PET Images#################################
##########################Averaged PET Images#################################
##########################Averaged PET Images#################################
AD_image_m <- apply(AD_image,2,mean)
AD1 <- as.vector(mask)
AD1[which(AD1!=0)] <- AD_image_m
AD1.arr <- array(AD1, dim=c(91, 109, 91))

MCItmp <- apply(MCI_image,2,mean) 
MCI1 <- as.vector(mask)
MCI1[which(MCI1!=0)] <- MCItmp
MCI1.arr <- array(MCI1, dim=c(91, 109, 91))

######Brodmann’s area 24/32/33
par(mfrow=c(2,2))
#X is right(+), left(-)
ImgX <- AD1.arr[47,,]
MyHeatMapADNI(ImgX)
ImgX <- AD1.arr[43,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- AD1.arr[,30,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- AD1.arr[,,45]
MyHeatMapADNI(ImgZ)
#X is right(+), left(-)
ImgX <- MCI1.arr[47,,]
MyHeatMapADNI(ImgX)
ImgX <- MCI1.arr[43,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- MCI1.arr[,30,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- MCI1.arr[,,45]
MyHeatMapADNI(ImgZ)

temp1 <- X_input[1:53,]
temp1 <- temp1[,colnames(temp1)%in%c(4011,4012)]
apply(temp1,2,mean)
temp2 <- X_input[-1:-53,]
temp2 <- temp2[,colnames(temp2)%in%c(4011,4012)]
apply(temp2,2,mean)

######Brodmann’s area 23/31
par(mfrow=c(2,2))
#X is right(+), left(-)
ImgX <- AD1.arr[48,,]
MyHeatMapADNI(ImgX)
ImgX <- AD1.arr[42,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- AD1.arr[,75,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- AD1.arr[,,50]
MyHeatMapADNI(ImgZ)
#X is right(+), left(-)
ImgX <- MCI1.arr[48,,]
MyHeatMapADNI(ImgX)
ImgX <- MCI1.arr[42,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- MCI1.arr[,75,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- MCI1.arr[,,50]
MyHeatMapADNI(ImgZ)

######Brodmann’s area 35
par(mfrow=c(2,2))
#X is right(+), left(-)
ImgX <- AD1.arr[40,,]
MyHeatMapADNI(ImgX)
ImgX <- AD1.arr[50,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- AD1.arr[,55,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- AD1.arr[,,60]
MyHeatMapADNI(ImgZ)
#X is right(+), left(-)
ImgX <- MCI1.arr[40,,]
MyHeatMapADNI(ImgX)
ImgX <- MCI1.arr[50,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- MCI1.arr[,55,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- MCI1.arr[,,60]
MyHeatMapADNI(ImgZ)

temp1 <- X_input[1:53,]
temp1 <- temp1[,colnames(temp1)%in%c(4101,4102)]
apply(temp1,2,mean)
temp2 <- X_input[-1:-53,]
temp2 <- temp2[,colnames(temp2)%in%c(4101,4102)]
apply(temp2,2,mean)

######Brodmann’s area 39
par(mfrow=c(2,2))
#X is right(+), left(-)
ImgX <- AD1.arr[25,,]
MyHeatMapADNI(ImgX)
ImgX <- AD1.arr[65,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- AD1.arr[,60,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- AD1.arr[,,26]
MyHeatMapADNI(ImgZ)
#X is right(+), left(-)
ImgX <- MCI1.arr[25,,]
MyHeatMapADNI(ImgX)
ImgX <- MCI1.arr[65,,]
MyHeatMapADNI(ImgX)
#Y is front(+) back(-)
ImgY <- MCI1.arr[,60,]
MyHeatMapADNI(ImgY)
#Z is top(+) down(-)
ImgZ <- MCI1.arr[,,26]
MyHeatMapADNI(ImgZ)

##########################Correlation matrices#################################
##########################Correlation matrices#################################
##########################Correlation matrices#################################
corRegionM_AD <- cor(X_input[1:53,])
corRegionM_AD[which(abs(corRegionM_AD)<0.8, arr.ind=TRUE)] <- 0
library(gplots)
heatmap.2(corRegionM_AD, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

corRegionM_MCI <- cor(X_input[-1:-53,])
corRegionM_MCI[which(abs(corRegionM_MCI)<0.8, arr.ind=TRUE)] <- 0
library(gplots)
heatmap.2(corRegionM_MCI, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

##########################Glass brain #################################
##########################Glass brain #################################
##########################Glass brain #################################
library(ggseg3d)
library(dplyr)
library(tidyr)

data(aseg_3d)
data(dk_3d)
comb_3d <- rbind(aseg_3d,dk_3d)

someData = comb_3d %>% 
  unnest(cols = ggseg_3d) %>% 
  select(label) %>% 
  filter(label=="rh_caudalanteriorcingulate"|label=="rh_posteriorcingulate"|label=="rh_rostralanteriorcingulate"|
           label=="rh_isthmuscingulate"|label=="Right-Hippocampus"|label=="rh_inferiorparietal"|label=="Right-Putamen") %>% 
  na.omit() %>% 
  mutate(p = seq(1, nrow(.)))

ggseg3d(.data = someData, 
        atlas = comb_3d,
        colour = "p", text = "p", 
        na.alpha= .2) %>% 
  add_glassbrain() %>% 
  remove_axes() 

##########################network#################################
##########################network#################################
##########################network#################################
networkRegionM_AD <- regionM[1:53,]
colnames(networkRegionM_AD) <- region_code[,3]
networkRegionM_AD <- round(cov(networkRegionM_AD),3)
network_AD <- solve(networkRegionM_AD)
network_AD[which(abs(corRegionM_AD)<0.8, arr.ind=TRUE)] <- 0
write.csv(network_AD,"network_AD.csv")
network_AD1 <- network_AD
network_AD1[which(network_AD1!=0)] <- 1
temp <- network_AD1[rownames(network_AD1)%in%c("Hippocampus_L","Hippocampus_R",
                                               "Cingulum_Mid_L","Cingulum_Mid_R","Cingulum_Post_L","Cingulum_Post_R",
                                               "Angular_L","Angular_R","Putamen_L","Putamen_R","ParaHippocampal_R",
                                               "Amygdala_R","Temporal_Mid_L"),
                    colnames(network_AD1)%in%c("Hippocampus_L","Hippocampus_R",
                                                "Cingulum_Mid_L","Cingulum_Mid_R","Cingulum_Post_L","Cingulum_Post_R",
                                                "Angular_L","Angular_R","Putamen_L","Putamen_R","ParaHippocampal_R",
                                                "Amygdala_R","Temporal_Mid_L")]
temp <- graph_from_adjacency_matrix(temp , mode='undirected', diag=F )
plot(temp, layout=layout.sphere, main="sphere")
plot(temp, layout=layout.circle, main="circle")
plot(temp, layout=layout.random, main="random")
plot(temp, layout=layout.fruchterman.reingold, main="Network of ROIs in AD group")
corRegionM_AD[which(abs(corRegionM_AD)<0.8, arr.ind=TRUE)] <- 0
library(gplots)
heatmap.2(corRegionM_AD, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

networkRegionM_MCI <- regionM[54:184,]
colnames(networkRegionM_MCI) <- region_code[,3]
networkRegionM_MCI <- round(cov(networkRegionM_MCI),3)
network_MCI <- solve(networkRegionM_MCI)
network_MCI[which(abs(corRegionM_MCI)<0.8, arr.ind=TRUE)] <- 0
write.csv(network_MCI,"network_MCI.csv")
network_MCI1 <- network_MCI
network_MCI1[which(network_MCI1!=0)] <- 1
temp <- network_MCI1[rownames(network_MCI1)%in%c("Hippocampus_L","Hippocampus_R",
                                                 "Cingulum_Mid_L","Cingulum_Mid_R","Cingulum_Post_L","Cingulum_Post_R",
                                                 "Angular_L","Angular_R","Putamen_L","Putamen_R","Precentral_L", "Precentral_R", "Frontal_Sup_L", "Frontal_Sup_R", "Postcentral_R",
                                                 "Frontal_Sup_R", "Frontal_Mid_R", "Frontal_Inf_Oper_R", "Cingulum_Ant_R",
                                                 "ParaHippocampal_L", "ParaHippocampal_R", "Amygdala_L", "Amygdala_R",
                                                 "Occipital_Mid_L", "Parietal_Inf_L", "SupraMarginal_L", "Temporal_Mid_L",
                                                 "Frontal_Mid_R", "Occipital_Sup_R", "Occipital_Mid_L", "Occipital_Mid_R", "Parietal_Inf_R",
                                                 "SupraMarginal_R", "Precuneus_R", "Temporal_Sup_R", "Temporal_Mid_R","Temporal_Inf_R"),
                    colnames(network_MCI1)%in%c("Hippocampus_L","Hippocampus_R",
                                               "Cingulum_Mid_L","Cingulum_Mid_R","Cingulum_Post_L","Cingulum_Post_R",
                                               "Angular_L","Angular_R","Putamen_L","Putamen_R","Precentral_L", "Precentral_R", "Frontal_Sup_L", "Frontal_Sup_R", "Postcentral_R",
                                               "Frontal_Sup_R", "Frontal_Mid_R", "Frontal_Inf_Oper_R", "Cingulum_Ant_R",
                                               "ParaHippocampal_L", "ParaHippocampal_R", "Amygdala_L", "Amygdala_R",
                                               "Occipital_Mid_L", "Parietal_Inf_L", "SupraMarginal_L", "Temporal_Mid_L",
                                               "Frontal_Mid_R", "Occipital_Sup_R", "Occipital_Mid_L", "Occipital_Mid_R", "Parietal_Inf_R",
                                               "SupraMarginal_R", "Precuneus_R", "Temporal_Sup_R", "Temporal_Mid_R","Temporal_Inf_R")]
temp <- graph_from_adjacency_matrix(temp , mode='undirected', diag=F )
temp <- network_MCI1
plot(temp, layout=layout.sphere, main="sphere")
plot(temp, layout=layout.circle, main="circle")
plot(temp, layout=layout.random, main="random")
plot(temp, layout=layout.fruchterman.reingold, main="Network of ROIs in AD group")


library(igraph)
temp <- pre_network1
temp <- graph_from_adjacency_matrix(temp , mode='undirected', diag=F )
#gg7 <- graph.edgelist(temp,directed=F) 
temp <- plot(temp,layout=layout.circle, vertex.size=2, vertex.label.dist=0.5,
             vertex.label.cex=.8,rescale=FALSE, xlim=c(-1,1), ylim=c(-1,1))

#####Circular network
library(edgebundleR)
# Create the graph object
temp <- pre_network1_1
g <- graph.data.frame(temp, directed=F)
edgebundle( g )

###################heatmap for network structure################
my_colors<- colorRampPalette(c("red", "green"))             # Manual color range
heatmap(pre_network1, col = my_colors(100))
my_colors<- colorRampPalette(c("red", "green"))             # Manual color range
heatmap(pre_network1, col = my_colors(2))

library(gplots)
heatmap.2(pre_network1, dendrogram='none', col=c("white", "red"), trace="none", key=FALSE)

pre_network2_HDnetGLM <- pre_network1
#Putamen_R
#Angular_R
#Cingulum_Post_R
#Paracentral_Lobule_L
#Hippocampus_L
#Cingulum_Mid_L
#Cingulum_Post_L
#Cingulum_Ant_L

pre_network2_HDnetGLM[which(rownames(pre_network2_HDnetGLM)=="Cingulum_Ant_L"),
                      which(colnames(pre_network2_HDnetGLM)=="Cingulum_Post_L") ] <- -1

pre_network2_HDnetGLM[which(rownames(pre_network2_HDnetGLM)%in%c("Hippocampus_L","Cingulum_Mid_L","Cingulum_Post_L","Cingulum_Ant_L","Putamen_R","Angular_R","Cingulum_Post_R","Paracentral_Lobule_L")),
                      which(colnames(pre_network2_HDnetGLM)%in%c("Hippocampus_L","Cingulum_Mid_L","Cingulum_Post_L","Cingulum_Ant_L","Putamen_R","Angular_R","Cingulum_Post_R","Paracentral_Lobule_L")) ]

heatmap.2(pre_network2_HDnetGLM, dendrogram='none', col=c("green","white", "red"), trace="none", key=FALSE)

pre_network2_Lasso <- pre_network1 
#Cingulum_Post_L
#Cingulum_Post_R
#Postcentral_L
#Angular_R
#Putamen_R
#Cerebelum_4_5_L

pre_network2_Lasso[which(rownames(pre_network2_Lasso)=="Cerebelum_4_5_L"),
                      which(colnames(pre_network2_Lasso)=="Cerebelum_4_5_L") ] <- -1

pre_network2_Lasso[which(rownames(pre_network2_Lasso)%in%c("Cingulum_Post_L", "Cingulum_Post_R", "Postcentral_L", "Angular_R", "Putamen_R", "Cerebelum_4_5_L")),
                      which(colnames(pre_network2_Lasso)%in%c("Cingulum_Post_L", "Cingulum_Post_R", "Postcentral_L", "Angular_R", "Putamen_R", "Cerebelum_4_5_L")) ]
heatmap.2(pre_network2_Lasso, dendrogram='none', col=c("green","white", "red"), trace="none", key=FALSE)
