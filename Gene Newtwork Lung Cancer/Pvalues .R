
################################################################
################################################################
#### Bootstrap-based permutation test on a selected marker "vox_168602"

#X_train[,All_set]

ISresample.12M <- NULL
ISresample.13M <- NULL
ISresample.23M <- NULL


for(subsample_id in 1:20000){


set.seed(subsample_id*4)

AD_subsample_id <- sample(c(1:40), 40, replace=TRUE)
MCI_subsample_id <- sample(c(41:144), 104, replace=TRUE)
NC_subsample_id <- sample(c(145:200), 56, replace=TRUE)


write.table(AD_subsample_id, file=paste("R", r, "_sample/AD_subsample_id", subsample_id, "_09282019.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
write.table(MCI_subsample_id, file=paste("R", r, "_sample/MCI_subsample_id", subsample_id, "_09282019.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)
write.table(NC_subsample_id, file=paste("R", r, "_sample/NC_subsample_id", subsample_id, "_09282019.txt", sep=""), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)


y_subsample <- c(All_y[AD_subsample_id], All_y[MCI_subsample_id], All_y[NC_subsample_id])
X_subsample <- X_train[c(AD_subsample_id, MCI_subsample_id, NC_subsample_id), All_set]
Z_subsample <- All_cova[c(AD_subsample_id, MCI_subsample_id, NC_subsample_id),2:3]

All_y <- y_train #y_subsample
X_All <- X_subsample
Z0 <- Z_subsample

X_All <- apply(X_All, 2, function(x) (x-mean(x))/sd(x))
Z0 <- apply(Z0, 2, function(x) (x-mean(x))/sd(x))

IS.M <- matrix(0, dim(X_All)[2], 3)

#mean12.Z <- apply(Z0[which(All_y==1),], 2, mean) - apply(Z0[which(All_y==2),], 2, mean)
#mean13.Z <- apply(Z0[which(All_y==1),], 2, mean) - apply(Z0[which(All_y==3),], 2, mean) 
#mean23.Z <- apply(Z0[which(All_y==2),], 2, mean) - apply(Z0[which(All_y==3),], 2, mean) 

#ZOmega <- qr.solve(cov(Z0))

#IS.M[,1] <- IS.M[,1] + ZOmega%*%mean12.Z
#IS.M[,2] <- IS.M[,2] + ZOmega%*%mean13.Z
#IS.M[,3] <- IS.M[,3] + ZOmega%*%mean23.Z

##

mean12.X <- apply(X_All[which(All_y==1),], 2, mean) - apply(X_All[which(All_y==2),], 2, mean)
mean13.X <- apply(X_All[which(All_y==1),], 2, mean) - apply(X_All[which(All_y==3),], 2, mean) 
mean23.X <- apply(X_All[which(All_y==2),], 2, mean) - apply(X_All[which(All_y==3),], 2, mean) 


Omega.r <- my.inv(cov(X_All))
       
IS.M[,1] <- IS.M[,1] + Omega.r%*%mean12.X
IS.M[,2] <- IS.M[,2] + Omega.r%*%mean13.X
IS.M[,3] <- IS.M[,3] + Omega.r%*%mean23.X


##

row.names(IS.M) <- All_vox

ISresample.12M <- cbind(ISresample.12M, IS.M[,1])
ISresample.13M <- cbind(ISresample.13M, IS.M[,2])
ISresample.23M <- cbind(ISresample.23M, IS.M[,3])

}



percentile <- ecdf(abs(ISresample.12M)[which(All_vox=="vox_168602"),])
1-percentile(abs(ISmatrix)[which(All_vox=="vox_168602"),1])
[1] 0




##################################################
##################################################
######
###### Hotelling’s test on a selected network


###############
###############
### from mLDA results
###

node_r <- names(network.nodes_list[[1]])  ######### 1 ADMCI,2 ADNC,3 MCINC 


Omega_r <- network.Omega_list[[1]]		######### 1 ADMCI,2 ADNC,3 MCINC 


MI_node <- as.vector(MI_matrix[,1])		######### 1 ADMCI,2 ADNC,3 MCINC 




Sele_node <- rownames(ISmatrix)[which(as.vector(ISmatrix[,1])!=0)] ######### 1 ADMCI,2 ADNC,3 MCINC 

Sele_node_IS <- ISmatrix[which(as.vector(ISmatrix[,1])!=0),1] ######### 1 ADMCI,2 ADNC,3 MCINC 

Sele_id <- match(Sele_node, node_r)
MI_id <- match(MI_node, node_r)


Adj <- abs(Omega_r)
colnames(Adj) <- node_r
rownames(Adj) <- node_r


vis=rep(0, dim(Adj)[1])
id <- match(voxel.name, rownames(Adj))

tmp <- con_node(id, Adj, visited=vis, m=dim(Adj)[1])$connect

node_r_new <- rownames(Adj)[tmp]
Adj_new <- Adj[tmp, tmp]

length(node_r_new)   #### network size

tmp0 <- match(Sele_node, node_r_new)
tmp0 <- tmp0[which(!is.na(tmp0))]
Sele_node_new <- node_r_new[tmp0]   
Sele_node_new         #### Selected voxels in the network
length(Sele_node_new) #### Number of selected voxels in the network

Sele_marg_id <- match(Sele_node_new, colnames(All_image))

tmp1 <- match(Sele_node_new, node_r)

Omega_r_sele <- Omega_r[tmp1, tmp1]

ADMCI.meanDiff_sele.v <-  ADMCI.meanDiff_ALL.v[Sele_marg_id]  ######### 1 ADMCI,2 ADNC,3 MCINC 



library(DescTools)


XpairSel <- All_image[which(All_y!=3), Sele_marg_id]
ypair <- All_y[which(All_y!=3)]

HotellingsT2Test(XpairSel ~ ypair)



	Hotelling's two sample T2-test

data:  XpairSel by ypair
T.2 = 4.7323, df1 = 24, df2 = 100, p-value = 1.644e-08
alternative hypothesis: true location difference is not equal to c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)


