trackgroupout <- function(FullTableFull,pooledvcov,testvars,nontestvars,sigmai,
    betai,betaWFE,TrackGroup,YYlist,XXlist) {
    #' Algorithm to move into closest group
    #'
    #' Algorithm to move into closest group
    #'
    #' @param FullTableFull data
    #' @return New groupings
    #' @export
    #' @examples
    #'
    
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    FullTableFull$Vessel <- as.factor(FullTableFull$Vessel)
    FullTableFull$EQ <- as.numeric(FullTableFull$EQ)
    
    for (i in levels(FullTableFull$Vessel)) {
        
        FullTable <- subset(FullTableFull, Vessel == i)
        
        testdi <- list()
        
        for (k in levels(FullTableFull$Group)) {
            
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
            
            testdi[[k]] <- abs((tempdi - length(betaitemp))/
			    (sqrt(2 * length(betaitemp))))

        }
        
        index <- match(Reduce("min", testdi), testdi)
        
        if (nrow(subset(TrackGroup, TrackGroup$Group ==
		    TrackGroup$Group[TrackGroup$Vessel == i])) > 1) {
			
            TrackGroup$Group[TrackGroup$Vessel == i] <- index
			
        }
        
    }
    
    return(TrackGroup)
    
}
