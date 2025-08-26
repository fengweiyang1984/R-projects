setwd("C:/Users/Fengwei Yang/Desktop/YFW work/GRA/Disertation 3/Dataset/LogisticModelData")
load("ADNI_info.Rdata")
#ls() function in R Language is used to list the names of all the objects that are present in the working directory.
ls()
#[1] "ADNI_AD_bl"    "ADNI_AD_m12"   "ADNI_AD_m6"    "ADNI_MCI_bl"  
#[5] "ADNI_MCI_m12"  "ADNI_MCI_m6"   "ADNI_NORM_bl"  "ADNI_NORM_m12"
#[9] "ADNI_NORM_m6" 
dim(ADNI_AD_bl)
#[1] 73 13
dim(ADNI_AD_m12)
#[1] 73 13
ADNI_AD_bl[1:10,]
#    Subject.ID Research.Group DX.Group visit.Time NPISCORE FAQ.Total.Score
#4   003_S_1059        Patient       AD   Baseline        2              23
#13  005_S_0221        Patient       AD   Baseline       13              22

#    MMSE.Total.Score CDR.Total.Score GDS.Total.Score Subject.Sex Subject.Age
#4                 NA              NA              NA           F       84.72
#13                NA              NA              NA           M       67.57

#    Age.Qualifier Subject.WeightKg
#4               Y            76.66
#13              Y           112.94

ADNI_AD_m12[1:10,]

load("ADNI_PET.RData")
ls()
# [1] "ADNI_AD_bl"    "ADNI_AD_m12"   "ADNI_AD_m6"    "ADNI_MCI_bl"  
# [5] "ADNI_MCI_m12"  "ADNI_MCI_m6"   "ADNI_NORM_bl"  "ADNI_NORM_m12"
# [9] "ADNI_NORM_m6"  "coord"         "dat"           "mask"         
#[13] "region_code"  

dim(coord)
#[1] 902629      3

names(dat)
#[1] "AD"   "MCI"  "NORM"

names(dat$AD)
#[1] "bl"  "m6"  "m12" "m24"
dim((dat$AD)$bl)
#[1]     54 185405




###############################################
###############################################
#### 03192019
####
####

load("ADNI_PET.RData")
ls()
#[1] "coord"       "dat"         "mask"        "region_code"

dim(mask)
#Loading required package: oro.nifti
#oro.nifti 0.9.1
#[1]  91 109  91

dim(coord)
#[1] 902629      3

coord[1:10, ]
#     x    y   z
#1  -88 -124 -70
#2  -86 -124 -70
#3  -84 -124 -70
#4  -82 -124 -70
#5  -80 -124 -70
#6  -78 -124 -70
#7  -76 -124 -70
#8  -74 -124 -70
#9  -72 -124 -70
#10 -70 -124 -70

91*91*109
#[1] 902629

ls(dat)
#[1] "AD"   "MCI"  "NORM"

ls(dat$AD)
#[1] "bl"  "m12" "m24" "m6" 

dim((dat$AD)$bl)
#[1]     54 185405

length(which(as.vector(mask)!=0))
#[1] 185405

#####
dim((dat$AD)$bl)
#[1]     54 185405

dim((dat$MCI)$bl)
#[1]    131 185405
 
dim((dat$NORM)$bl)
#[1]     72 185405




AD_b <- (dat$AD)$bl
dim(AD_b)
#[1]     54 185405
write.table(AD_b, file="AD_bl.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)


MCI_b <- (dat$MCI)$bl
dim(MCI_b)
#[1]    131 185405
write.table(MCI_b, file="MCI_bl.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)


NORM_b <- (dat$NORM)$bl
dim(NORM_b)
#[1]     72 185405
write.table(NORM_b, file="NORM_bl.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)



dim(ADNI_AD_bl)
#[1] 73 13
as.vector(ADNI_AD_bl[,1])
as.vector(rownames(AD_b))
match(as.vector(rownames(AD_b)), as.vector(ADNI_AD_bl[,1]))


dim(ADNI_MCI_bl)
#[1] 172  13
as.vector(ADNI_MCI_bl[,1])
as.vector(rownames(MCI_b))
match(as.vector(rownames(MCI_b)), as.vector(ADNI_MCI_bl[,1]))


dim(ADNI_NORM_bl)
#[1] 81 13
as.vector(ADNI_NORM_bl[,1])
as.vector(rownames(NORM_b))
match(as.vector(rownames(NORM_b)), as.vector(ADNI_NORM_bl[,1]))


###

AD_mid <- match(as.vector(rownames(AD_b)), as.vector(ADNI_AD_bl[,1]))
AD_image <- AD_b[which(!is.na(AD_mid)), ]
AD_info <- ADNI_AD_bl[AD_mid[which(!is.na(AD_mid))],]

write.table(AD_image, file="AD_bl_image.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
write.table(AD_info, file="AD_bl_info.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

MCI_mid <- match(as.vector(rownames(MCI_b)), as.vector(ADNI_MCI_bl[,1]))

MCI_image <- MCI_b[which(!is.na(MCI_mid)), ]
MCI_info <- ADNI_MCI_bl[MCI_mid[which(!is.na(MCI_mid))],]

write.table(MCI_image, file="MCI_bl_image.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
write.table(MCI_info, file="MCI_bl_info.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

NORM_mid <- match(as.vector(rownames(NORM_b)), as.vector(ADNI_NORM_bl[,1]))

NORM_image <- NORM_b[which(!is.na(NORM_mid)), ]
NORM_info <- ADNI_NORM_bl[NORM_mid[which(!is.na(NORM_mid))],]

write.table(NORM_image, file="NORM_bl_image.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=TRUE)
write.table(NORM_info, file="NORM_bl_info.txt", sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)


###








dim(mask)
#[1]  91 109  91

ImgZ <- mask[,,45]
MyHeatMapADNI <- function(Beta){
  P <-nrow(Beta)
  Q <-ncol(Beta)
  Xlim <- c(0,2*(P+1))
  Ylim <- c(0,2*(Q+1))
  RGBColors <- col2rgb(colors()[1:length(colors())])
  HSVColors <- rgb2hsv( RGBColors[1,], RGBColors[2,], RGBColors[3,], maxColorValue=255)
  HueOrder <- order( HSVColors[1,], HSVColors[2,], HSVColors[3,] )
  uBeta <- unique(as.vector(Beta))
  ruBeta <- rank(uBeta)
  vect <- cbind(uBeta, ruBeta)
  l <- length(unique(ruBeta))
  plot(0, type="n", xlab="", ylab="", xlim=Xlim, ylim=Ylim, cex.lab=1.0, bty="n", axes=F)
 for (p in 1:P){
 for (q in 1:Q){
         k0 <- ruBeta[which(uBeta==Beta[p,q])]
 k <- round(k0*212/l)
 if(k==0){ k <-1}
 rect(2*(P-p+1)-1,2*(Q-q+1)-1, 2*(P-p+1)+1, 2*(Q-q+1)+1, col=colors()[HueOrder[k]], border=NA, lwd=0)
}
}
}

MyHeatMapADNI(ImgZ)

ImgY <- mask[,54,]
MyHeatMapADNI(ImgY)



ADtmp <- AD_b[1,] 
AD1 <- as.vector(mask)
AD1[which(AD1!=0)] <- ADtmp
AD1.arr <- array(AD1, dim=c(91, 109, 91))
 
dim(AD1.arr)
#[1]  91 109  91

ImgZ <- AD1.arr[,,45]
#> X11()
MyHeatMapADNI(ImgZ)

ImgY <- AD1.arr[,54,]
MyHeatMapADNI(ImgY)

NORMtmp <- NORM_b[1,] 
NORM1 <- as.vector(mask)
NORM1[which(NORM1!=0)] <- NORMtmp
NORM1.arr <- array(NORM1, dim=c(91, 109, 91))

dim(NORM1.arr)
#[1]  91 109  91

ImgZ <- NORM1.arr[,,45]
X11()
MyHeatMapADNI(ImgZ)

ImgY <- NORM1.arr[,54,]
MyHeatMapADNI(ImgY)

MCItmp <- MCI_b[1,] 
MCI1 <- as.vector(mask)
MCI1[which(MCI1!=0)] <- MCItmp
MCI1.arr <- array(MCI1, dim=c(91, 109, 91))

All_b <- rbind(AD_b, MCI_b, NORM_b)

All_image <- rbind(AD_image, MCI_image, NORM_image) 
All_info <- rbind(AD_info, MCI_info, NORM_info)
All_info$Subject.Sex.Num <- (All_info$Subject.Sex=="M")*1


mask_non0 <- as.vector(mask)[which(as.vector(mask)!=0)]
table(mask_non0)
#mask_non0
#2001 2002 2101 2102 2111 2112 2201 2202 2211 2212 2301 2302 2311 2312 2321 2322 2331 2332 2401 2402 2501 2502 2601 2602 2611 2612 
#3526 3381 3599 4056  963  997 4863 5104  888 1015 1038 1399 2529 2151 1690 1707  990 1331 2147 2371  280  289 2992 2134  719  856 
#2701 2702 3001 3002 4001 4002 4011 4012 4021 4022 4101 4102 4111 4112 4201 4202 5001 5002 5011 5012 5021 5022 5101 5102 5201 5202 
#852  745 1858 1770 1400 1313 1941 2203  463  335  932  946  978 1132  220  248 2258 1861 1526 1424 2095 2300 1366 1413 3270 2098 
#5301 5302 5401 5402 6001 6002 6101 6102 6201 6202 6211 6212 6221 6222 6301 6302 6401 6402 7001 7002 7011 7012 7021 7022 7101 7102 
#941  989 2310 2518 3892 3823 2065 2222 2447 1345 1256 1974 1173 1752 3528 3265 1349  836  962  994 1009 1064  293  280 1100 1057 
#8101 8102 8111 8112 8121 8122 8201 8202 8211 8212 8301 8302 9001 9002 9011 9012 9021 9022 9031 9032 9041 9042 9051 9052 9061 9062 
#225  249 2296 3141 1285 1338 4942 4409  755 1187 3200 3557 2603 2648 1894 2117  136  207 1125  861 1694 1795  585  534 1887 2308 
#9071 9072 9081 9082 9100 9110 9120 9130 9140 9150 9160 9170 
#869  809  144  159   53  228  665  371  194  243  174  112 

Tb_mask_non0 <- table(mask_non0)

regionM <- NULL
for (t in 1:length(Tb_mask_non0)){
	tmpM <- All_image[, which(mask_non0==names(Tb_mask_non0)[t])]
	tmp.v <- apply(tmpM, 1, mean)
	regionM <- cbind(regionM, tmp.v)
}

colnames(regionM) <- names(Tb_mask_non0)

corRegionM <- cor(regionM)
corRegionM0 <- corRegionM

corRegionM <- corRegionM0
corRegionM[which(abs(corRegionM)<0.8, arr.ind=TRUE)] <- 0


library(gplots)
#pdf("threshCorrMatrix08.pdf")
png("threshCorrMatrix08.png")
heatmap.2(corRegionM, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')
dev.off()


write.table(corRegionM, file="threshCorrMatrix08.txt", sep="\t", quote=FALSE, row.names=TRUE, col.names=TRUE)

#diag(corRegionM) <- 0
#regionIncl <- apply(corRegionM, 2, sum)
#regionIncl

All_cova <- All_info[,c(10, 11, 13, 14)]
All_y <- c(rep(1, dim(AD_image)[1]), rep(2, dim(MCI_image)[1]), rep(3, dim(NORM_image)[1])) 

r <- 1
corRegion.v0 <- corRegionM[r,]
corRegion.v <- which(corRegion.v0!=0)

ImageM.r <- NULL
for(j in 1:length(corRegion.v)){

regionID <- which(mask_non0==as.numeric(names(corRegion.v)[1]))
	ImageM.r <- cbind(ImageM.r, All_image[,regionID])
}

tmp1 <- cor(ImageM.r)
heatmap.2(tmp1, dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none')

#what's this???
region_code


################################################################
################################################################
mLDA.Kclass <-
function(X, y, X.new, Z=NULL, COV=NULL, COR=NULL, K=3, tau=20, alpha=0.5, nu=100, d=10){

    p <- dim(X)[2]
    n <- dim(X)[1]

    if(length(y)!=n){
        stop("X and Y contain different number of samples!!")
     }
    if(dim(X.new)[2]!=p){
        stop("The test data do not contain the same number of features as the training data!!")
     }
    
    if(is.null(COV)){
	COV <- cov(X)
    }
  
    if(is.null(COR)){
	COR <- cor(X)
    }


    ISmatrix <- matrix(0, p, (K*(K-1)/2))

	MImatrix <- matrix(0, tau, (K*(K-1)/2))
	conn_Array <- array(0, dim=c(tau, p, (K*(K-1)/2)))

    for(k1 in 2:K){
   	 for(k in 1:(k1-1)){
   	    tmp <- mLDA.pair(X, y, X.new, Z=Z, COV=COV, COR=COR, pair=c(k,k1), tau=tau, alpha=alpha, nu=nu, d=d)     
 		ISmatrix[as.numeric(names(tmp$iffcond)), ((k1-1)*(k1-2)/2+k)] <- tmp$iffcond
		MImatrix[,((k1-1)*(k1-2)/2+k)] <- tmp$MIset
		conn_Array[,,((k1-1)*(k1-2)/2+k)] <- tmp$connMatrix
	  }
    }

    nu.vec.tmp <- NULL
    for(k in 1:K){
         nu.vec.tmp <- c(nu.vec.tmp, ISmatrix[which(ISmatrix[,1]!=0), k])
    }

    nu.value <- sort(nu.vec.tmp, decreasing=TRUE)[nu]

    ISmatrix[which(ISmatrix<nu.value, arr.ind=TRUE)] <- 0
    Sset <- which(apply(ISmatrix, 1, function(x) sum(abs(x)))!=0)
    ISmatrix.o <- ISmatrix[Sset, ]

    Sset.tmp <- apply(ISmatrix, 1, function(x) sum(abs(x)))[Sset]
    names(Sset.tmp) <- Sset
    Sset.order <- sort(Sset.tmp, decreasing=TRUE)

    mean.cis.mat <- NULL

    for(k in 1:K){
       id <- which(y==k)
       X.k <- X[id, Sset]
       mean.cis.vec <- apply(X.k, 2, mean)
       mean.cis.mat <- cbind(mean.cis.mat, mean.cis.vec)
    }

    COV.cis.mat <- COV[Sset, Sset]
    Omega.cis.mat <-  my.inv(COV.cis.mat)

    Fisher.matrix <- NULL

    for(i in 1:dim(X.new)[1]){
       X.new.cis.i <- X.new[i, Sset]
       Fisher.vec <- NULL
       for(k in 1:K){
           Fisher.ik <- (X.new.cis.i-mean.cis.mat[,k]/2)%*%Omega.cis.mat%*%mean.cis.mat[,k]
           Fisher.vec <- c(Fisher.vec, Fisher.ik)
       }
       Fisher.matrix <- rbind(Fisher.matrix, Fisher.vec)
    }

    row.names(Fisher.matrix) <- c(1:(dim(X.new)[1]))

    predClass <- apply(Fisher.matrix, 1, which.max)


    result <- list(ISmatrix=ISmatrix.o, screenset=Sset, Fisher.matrix=Fisher.matrix, PredClass=predClass, MImatrix=MImatrix, conn_Array=conn_Array, Omega.cis.mat = Omega.cis.mat)

    result


}










mLDA.pair <-
function(X, y, X.new, Z=NULL, COV=NULL, COR=NULL, pair=c(1,2), tau=20, alpha=0.5, nu=100, d=10){
  
    p <- dim(X)[2]
    n <- dim(X)[1]

    if(length(y)!=n){
        stop("X and Y contain different number of samples!!")
     }
    if(dim(X.new)[2]!=p){
        stop("The test data do not contain the same number of features as the training data!!")
     }

    id1<- which(y==pair[1])
    if(length(id1)==0){
		stop(paste("There is no y entries labeled as ", pair[1], "!", sep=""))
	}

    id2<- which(y==pair[2])
    if(length(id2)==0){
		stop(paste("There is no y entries labeled as ", pair[2], "!", sep=""))
	}

	y12 <- c(rep(0, length(id1)), rep(1, length(id2)))
	X12 <- X[c(id1, id2),]
    
    if(is.null(COV)){
		COV <- cov(X)
    }
  
    if(is.null(COR)){
		COR <- cor(X)
    }


	COR[which(abs(COR)<alpha, arr.ind=TRUE)] <- 0


	meanDiff <- apply(as.matrix(X12[which(y12==0),]), 2, mean) - apply(as.matrix(X12[which(y12==1),]), 2, mean)
	meanDiffRank <- order(-abs(meanDiff))
	meanDiff_firstTau <- meanDiffRank[1:tau]

	meanAvg <- (apply(as.matrix(X12[which(y12==0),]), 2, mean) + apply(as.matrix(X12[which(y12==1),]), 2, mean))/2

	Omega <- matrix(0, p, p)

	vis=rep(0, p)
	id <- meanDiff_firstTau[1]
	scr <- meanDiff_firstTau
	conn_matrix	<- matrix(0, tau, p)
    rownames(conn_matrix) <- meanDiff_firstTau
    colnames(conn_matrix) <- c(1:p)

	while(!is.na(id)){
		tmp <- con_node(id, COR, visited=vis, m=d)
		con_set <- tmp$connect
		vis <- tmp$visited
		Omega[con_set, con_set] <- qr.solve(COV[con_set, con_set])

		tau_id <- which(meanDiff_firstTau==id)
		conn_matrix[tau_id, con_set] <- 1 

		id <- scr[is.na(match(scr, con_set))][1]
		scr <- scr[is.na(match(scr, con_set))]
	}

	screenID <- which(apply(Omega, 1, sum)!=0)
	Omega.s <- Omega[screenID, screenID]

	if(!is.null(Z)){
#		Z <- apply(Z, 2, function(x) (x-mean(x))/sd(x))
		ZmeanDiff <- apply(as.matrix(Z[which(y12==0),]), 2, mean) - apply(as.matrix(Z[which(y12==1),]), 2, mean)
		ZmeanAvg <- (apply(as.matrix(Z[which(y12==0),]), 2, mean) + apply(as.matrix(Z[which(y12==1),]), 2, mean))/2
		ZOmega <- qr.solve(cov(Z))
		pZ <- dim(Z)[2]
	}else{
		ZmeanDiff <- NULL
		ZmeanAvg <- NULL
		ZOmega <- NULL
		pZ <- 0
	}

	meanDiff.s <- meanDiff[screenID]
	iffcond.s <- as.vector(Omega.s%*%meanDiff.s)
	names(iffcond.s) <- screenID
	iffcond.s_order <- sort(abs(iffcond.s), decreasing=TRUE)

    if(length(iffcond.s_order)>nu){
		screenID2 <- as.numeric(names(iffcond.s_order)[1:nu])
	}else{
 		screenID2 <- as.numeric(names(iffcond.s_order))
    }

	meanAvg.s0 <- meanAvg[screenID2] 
	Omega.s20 <- Omega[screenID2, screenID2] 
	meanDiff.s20 <- meanDiff[screenID2]
	X.new.s0 <- X.new[,screenID2]


	meanAvg.s <- c(ZmeanDiff, meanAvg.s0)
    pXs2 <- length(screenID2)
	Omega.s2 <- matrix(0, (pZ+pXs2), (pZ+pXs2)) 
	Omega.s2[1:pZ, 1:pZ] <- ZOmega
	Omega.s2[(pZ+1):(pZ+pXs2), (pZ+1):(pZ+pXs2)] <- Omega[screenID2, screenID2] 
#	Omega.s2 <- Omega[screenID2, screenID2] 
	meanDiff.s2 <- c(ZmeanDiff, meanDiff[screenID2])
	X.new.s <- cbind(Z, X.new[,screenID2])


	if(length(meanAvg.s)==1){
		meanAvg.s <- as.vector(meanAvg.s)
		X.new.s <- as.matrix(X.new.s)
		meanAvg.s.m <- matrix(rep(meanAvg.s, dim(X.new.s)[1]), dim(X.new.s)[1], length(meanAvg.s), byrow=TRUE)
	}else{
		meanAvg.s.m <- matrix(rep(meanAvg.s, dim(X.new.s)[1]), dim(X.new.s)[1], length(meanAvg.s), byrow=TRUE)
	}

#	ZmeanAvg.m0 <- matrix(rep(ZmeanAvg, dim(Z)[1]), dim(Z)[1], dim(Z)[2], byrow=TRUE)
#	meanAvg.s.m <- cbind(ZmeanAvg.m0, meanAvg.s.m0)


	FisherDR <- (X.new.s-meanAvg.s.m)%*%Omega.s2%*%meanDiff.s2
 	PredClass <- (FisherDR>=0)

	result <- list(iffcond=iffcond.s_order, screenset=screenID2, FisherDR=FisherDR, PredClass=PredClass, MIset = meanDiff_firstTau, connMatrix = conn_matrix)

    result
}











con_node <-
function(vid, Sig, visited=rep(0, dim(Sig)[1]), depth=0, m=10, connect=NULL, path=NULL){
         
     if(dim(Sig)[1] != dim(Sig)[2]){
        stop("Sig is not a square matrix")
     }
     if(vid>dim(Sig)[1] | vid<=0){
        stop("Wrong vertex index")
     }
  
     connect <- c(connect, vid)
     path <- c(path, vid)
     visited[vid] <- 1

     flag <- 0

     tmp <- which(visited==0)  
	 ranks <- order(-abs(Sig[vid,tmp]))
     rest.ranks <- tmp[ranks]

     for(i in rest.ranks){
       if(Sig[vid,i]!=0 & visited[i]==0 & depth<=m){
	   depth <- depth+1
	   tmp <- con_node(i, Sig, visited=visited, depth=depth, connect=connect, path=path)

	   visited <- tmp$visited
       depth <- tmp$depth
       connect <- tmp$connect
       path <- tmp$path
           flag <-1
        } # end of if        
     } # end of for

     if(flag==0){stop}

     return(list(connect=connect, visited=visited, depth=depth, path=path))
}











my.inv <-
function(X, eps=1e-12){
	eig.X <- eigen(X, symmetric=TRUE)
	P <- eig.X[[2]]
	lambda <- eig.X[[1]]
	ind <- lambda > eps
	lambda[ind] <- 1/lambda[ind]
	lambda[!ind] <- 0
	ans <- P%*%diag(lambda, nrow=length(lambda))%*%t(P)
	return(ans)
}









