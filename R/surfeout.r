surfeout <- function(FullTableFull,pooledvcov,testvars,nontestvars,YYlist,
    XXlist) {
    #' FE estimates for construction of dispersion test
    #'
    #' FE estimates for use in constructing dispersion test, gets called as part
    #' of constructing the test in kmeans_inner
    #'
    #' @param FullTableFull data
    #' @return SUR parameter estimates
    #' @export
    #' @examples
    #'
    
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    FullTableFull$Vessel <- as.factor(FullTableFull$Vessel)
    FullTableFull$EQ <- as.numeric(FullTableFull$EQ)
    
    SURbetaFE <- list()
    
    for (k in levels(FullTableFull$Group)) {
        
        FullTable2 <- subset(FullTableFull, Group == k)
        FullTable2$Vessel <- droplevels(FullTable2$Vessel)
        
        FEpart1 <- list()
        FEpart2 <- list()
        
        for (i in levels(FullTable2$Vessel)) {
            
            FullTable <- subset(FullTable2, Vessel == i)
            
            XX <- as.matrix(FullTable[unique(unlist(XXlist))])
            
            YY <- as.matrix(rowSums(FullTable[unique(unlist(YYlist))]))
            
            pooledBigSigma <- kronecker(solve(pooledvcov[[k]]),
			    diag((dim(YY)[1]/4)))
            cholPBS <- chol(pooledBigSigma)
			
			#note here that t(cholPBS %*% XX) %*% (cholPBS %*% XX) and
			    #t(XX) %*% (pooledBigSigma) %*% (XX) are numerically identical
				#(e.g. try all.equal()).
				#If applying FWL the order matters however, and need to
				#premultiply first.
            PcholXX <- cholPBS %*% XX
            PcholYY <- cholPBS %*% YY
            
			#see remarks 3 and 4 for form of partitioned unbalanced dispersion 
			    #tests in Pesaran & Yamagata 2008
            PcholXXtest <- PcholXX[, testvars]
            PcholXXnontest <- PcholXX[, nontestvars]
            PcholPartM <- diag(dim(YY)[1]) - PcholXXnontest %*%
			    solve(t(PcholXXnontest) %*% PcholXXnontest) %*%
				t(PcholXXnontest)
            
            FEpart1[[i]] <- (t(PcholXXtest) %*% PcholPartM %*% PcholXXtest)
            FEpart2[[i]] <- (t(PcholXXtest) %*% PcholPartM %*% PcholYY)
            
        }
        
        SURbetaFE[[k]] <- solve(Reduce("+", FEpart1)) %*% Reduce("+", FEpart2)
        
    }
    
    return(SURbetaFE)
    
}
