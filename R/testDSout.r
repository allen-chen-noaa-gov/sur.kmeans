testDSout <- function(FullTableFull,pooledvcov,testvars,nontestvars,sigmai,
    betai,betaWFE,Vessel,YYlist,XXlist) {
    #' Parameter homogeneity test
    #'
    #' Parameter homogeneity test
    #'
    #' @param FullTableFull data
    #' @return Parameter homogeneity test statistic
    #' @export
    #' @examples
    #'
    
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    FullTableFull$Vessel <- as.factor(FullTableFull$Vessel)
    FullTableFull$EQ <- as.numeric(FullTableFull$EQ)
    
    testDS <- list()
    
    for (k in levels(FullTableFull$Group)) {
        
        FullTable2 <- subset(FullTableFull, Group == k)
        FullTable2$Vessel <- droplevels(FullTable2$Vessel)
        
        testdi <- list()
        
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
            
            tempsigma <- sigmai[[i]]
            
            betaitemp <- betai[[i]]
            
            tempdi <- t(betaitemp - betaWFE[[k]]) %*% ((t(PcholXXtest) %*%
                PcholPartM %*% PcholXXtest)/as.numeric(tempsigma)) %*%
				(betaitemp - betaWFE[[k]])
            
            testdi[[i]] <- (tempdi - length(betaitemp))/
			    (sqrt(2 * length(betaitemp)))
            
        }
        
        testDS[[k]] <- (Reduce("+", testdi)) * (1/sqrt(length(Vessel)))
        
    }
    
    return(testDS)
    
}
