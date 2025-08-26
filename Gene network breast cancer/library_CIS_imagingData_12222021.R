################################################################
################################################################


mLDA.Kclass <-
function(X, y, X.new, Z=NULL, Z.new=NULL, COV=NULL, COR=NULL, K=3, tau=20, alpha=0.5, nu=100, d=10){

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


    marker.name.vec <- colnames(X)
    if(is.null(marker.name.vec)){
	marker.name.vec <- c(1:p) #paste("V", c(1:p), sep="")
	colnames(X) <- marker.name.vec
    }
    colnames(COV) <- marker.name.vec
    rownames(COV) <- marker.name.vec
    colnames(COR) <- marker.name.vec
    rownames(COR) <- marker.name.vec
  

	
    ISmatrix <- matrix(0, p, (K*(K-1)/2))
    rownames(ISmatrix) <- marker.name.vec

    MImatrix.name <- matrix(0, tau, (K*(K-1)/2))
    MImatrix.cord <- matrix(0, tau, (K*(K-1)/2))

    conn_Array <- array(0, dim=c(tau, p, (K*(K-1)/2)))
    dimnames(conn_Array)[[2]] <- marker.name.vec

    network.nodes_list <- NULL
    network.Omega_list <- NULL 
    
    local.network_LIST <- NULL

    for(k1 in 2:K){
   	 for(k in 1:(k1-1)){
   	    tmp <- mLDA.pair(X, y, X.new, Z=Z, Z.new=Z.new, COV=COV, COR=COR, pair=c(k,k1), tau=tau, alpha=alpha, nu=nu, d=d)     
 		ISmatrix[as.numeric(names(tmp$iffcond)), ((k1-1)*(k1-2)/2+k)] <- tmp$iffcond
		MImatrix.cord[,((k1-1)*(k1-2)/2+k)] <- tmp$MIset
		MImatrix.name[,((k1-1)*(k1-2)/2+k)] <- names(tmp$MIset)
		conn_Array[,,((k1-1)*(k1-2)/2+k)] <- tmp$connMatrix
		network.nodes_list[[((k1-1)*(k1-2)/2+k)]] <- tmp$network.nodes
		network.Omega_list[[((k1-1)*(k1-2)/2+k)]] <- tmp$network.Omega
		
		local.network_LIST[[((k1-1)*(k1-2)/2+k)]] <- tmp$local.network_list
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

    #### this is not really useful 08212020
#    Sset.tmp <- apply(ISmatrix, 1, function(x) sum(abs(x)))[Sset]
#    names(Sset.tmp) <- Sset
#    Sset.order <- sort(Sset.tmp, decreasing=TRUE)
    ########################################

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


    result <- list(network.nodes_list=network.nodes_list, network.Omega_list=network.Omega_list, local.network_LIST=local.network_LIST, ISmatrix=ISmatrix.o, screenset=Sset, Fisher.matrix=Fisher.matrix, PredClass=predClass, MImatrix.cord=MImatrix.cord, MImatrix.name=MImatrix.name, conn_Array=conn_Array)

    result


}

##### ISmatrix lists all "iffcond"s for only selected features in network.nodes_list (in an descent order).








mLDA.pair <-
function(X, y, X.new, Z=NULL, Z.new=NULL, COV=NULL, COR=NULL, pair=c(1,2), tau=20, alpha=0.5, nu=100, d=2, nb=10){
  
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


    marker.name.vec <- colnames(X)
    if(is.null(marker.name.vec)){
	marker.name.vec <- c(1:p) #paste("V", c(1:p), sep="")
	colnames(X) <- marker.name.vec
    }
    colnames(COV) <- marker.name.vec
    rownames(COV) <- marker.name.vec
    colnames(COR) <- marker.name.vec
    rownames(COR) <- marker.name.vec
#    colnames(X12) <- marker.name.vec

	COR[which(abs(COR)<alpha, arr.ind=TRUE)] <- 0


	meanDiff <- (apply(as.matrix(X12[which(y12==0),]), 2, mean) - apply(as.matrix(X12[which(y12==1),]), 2, mean))/apply(as.matrix(X12), 2, sd)
	#names(meanDiff) <- marker.name.vec
	meanDiffRank <- order(-abs(meanDiff))
	meanDiff_firstTau <- meanDiffRank[1:tau]
	names(meanDiff_firstTau) <- marker.name.vec[meanDiff_firstTau]

	meanAvg <- (apply(as.matrix(X12[which(y12==0),]), 2, mean) + apply(as.matrix(X12[which(y12==1),]), 2, mean))/2

	Omega <- matrix(0, p, p)
	colnames(Omega) <- marker.name.vec
	rownames(Omega) <- marker.name.vec

	vis=rep(0, p)
	id <- meanDiff_firstTau[1]
	scr <- meanDiff_firstTau
	conn_matrix	<- matrix(0, tau, p)
	rownames(conn_matrix) <- meanDiff_firstTau
	colnames(conn_matrix) <- marker.name.vec
	
	local.network_list <- list() 

	while(!is.na(id)){
		#tmp <- con_node(id, COR, visited=vis, m=d)
		tmp <- con_node_dn(id, COR, visited=vis, depth.max=d, neighbor.max=nb)
		con_set <- tmp$connect
		names(con_set) <- marker.name.vec[con_set]
		vis <- tmp$visited
		Omega[con_set, con_set] <- my.inv(COV[con_set, con_set])
		
		local.network_list <- append(local.network_list, list(con_set))

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
	iffcond.s_order_id <- order(abs(iffcond.s), decreasing=TRUE)
	iffcond.s_order <- iffcond.s[iffcond.s_order_id]

    if(length(iffcond.s_order)>nu){
		screenID2 <- as.numeric(names(iffcond.s_order)[1:nu])
	}else{
 		screenID2 <- as.numeric(names(iffcond.s_order))
    }

	meanAvg.s0 <- meanAvg[screenID2] 
	Omega.s20 <- Omega[screenID2, screenID2] 
	meanDiff.s20 <- meanDiff[screenID2]
	X.new.s0 <- X.new[,screenID2]

	if(!is.null(Z)){
		meanAvg.s <- c(ZmeanDiff, meanAvg.s0)
    	pXs2 <- length(screenID2)
		Omega.s2 <- matrix(0, (pZ+pXs2), (pZ+pXs2)) 
		Omega.s2[1:pZ, 1:pZ] <- ZOmega
		Omega.s2[(pZ+1):(pZ+pXs2), (pZ+1):(pZ+pXs2)] <- Omega[screenID2, screenID2] 
#	Omega.s2 <- Omega[screenID2, screenID2] 
		meanDiff.s2 <- c(ZmeanDiff, meanDiff[screenID2])
		X.new.s <- cbind(Z.new, X.new[,screenID2])
	}else{
		meanAvg.s <- meanAvg.s0
		Omega.s2 <- Omega[screenID2, screenID2] 
		meanDiff.s2 <- meanDiff[screenID2]
		X.new.s <- X.new[,screenID2]
	}

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

	result <- list(network.nodes=screenID, network.Omega=Omega.s, local.network_list=local.network_list, iffcond=iffcond.s_order, screenset=screenID2, FisherDR=FisherDR, PredClass=PredClass, MIset = meanDiff_firstTau, connMatrix = conn_matrix)

    result
}


##### iffcond lists iffcond for the first nu selected features in network.nodes (in an descent order). screenset lists only the first nu selected features.






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






####################
####################
### 08/06/2021
###




#==============================  final version of neighbor+layer selection
con_node_dn <-
  function(vid, Sig, visited=rep(0, dim(Sig)[1]), depth=0, depth.max=2, neighbor.max=10, connect=NULL, path=NULL, depth.path=NULL){
    
    if(dim(Sig)[1] != dim(Sig)[2]){ stop("Sig is not a square matrix") }
    if(vid>dim(Sig)[1] | vid<=0){ stop("Wrong vertex index")}

    neighbor <- 0
    
    if(neighbor<=neighbor.max & depth<=depth.max){    
    	connect <- c(connect, vid)
    	path <- c(path, vid)
    	visited[vid] <- 1
     }
    
    
    flag <- 0
    
    tmp <- which(visited==0)  
    ranks <- order(-abs(Sig[vid,tmp]))
    rest.ranks <- tmp[ranks]
    
    depth0 <- depth
    neighbor0 <- neighbor

    
    for(i in rest.ranks){
          depth1 <- depth0
             
      if(Sig[vid,i]!=0 & visited[i]==0 & neighbor<neighbor.max & depth1<depth.max){
      
      
      
        	neighbor <- neighbor +1

         
                depth1 <- depth1+1
         
        	tmp <- con_node_dn(i, Sig, visited=visited, depth=depth1, depth.max=depth.max, neighbor.max=neighbor.max, connect=connect, path=path, depth.path=depth.path)
        
        	visited <- tmp$visited
        	depth <- tmp$depth
         	connect <- tmp$connect
        	path <- tmp$path
        	
        	depth.path <- tmp$depth.path
        	depth.path <- c(depth1, depth.path)
        	
        	flag <-1

     }  # end of  if(Sig[vid,i]!=0 & visited[i]==0 & neighbor<=neighbor.max)
    } # end of for
    

    
   if(flag==0){stop}
    
    return(list(connect=connect, visited=visited, depth=depth, path=path, depth.path=depth.path))
  }





