kmeans_inner <- function(nn,FullTableFull,ntreat,Vessel,testvars,nontestvars,
    YYlist,XXlist) {
    #' K-means algorithm inner
    #'
    #' Keeps regrouping until convergence (original) per seed
    #'
    #' @param FullTableFull data
    #' @return Converged groupings and associated statistical tests
    #' @export
    #' @examples
    #'
    
    j <- nn
    
    set.seed(j)
    print(paste("Seed", j))
    
    Group <- as.factor(sample(1:ntreat, length(Vessel), replace = TRUE))
	
    while (min(table(Group)) < 1 || length(levels(Group)) != ntreat) {
	    Group <- as.factor(sample(1:ntreat, length(Vessel), replace = TRUE))
    }
    
    grouping <- cbind(Vessel, Group)
    
    Initgroups <- as.data.frame(grouping)
    Initgroups$Group <- as.numeric(as.character(Initgroups$Group))
    Initgroups2 <- Initgroups
    Difference <- rep(1, dim(Initgroups)[1])
    
    counting <- 1
    
    while (sum(abs(Difference)) > 0) {
        
        FullTableFull$Group <- NULL
        FullTableFull <- merge(FullTableFull, Initgroups, by = c("Vessel"))
		FullTableFull <- FullTableFull[order(FullTableFull$obsnum),]
        FullTableFull$Group <- as.factor(FullTableFull$Group)
        
        TrackGroup <- Initgroups
        
        betaFE <- indEQFEout(FullTableFull,YYlist,XXlist)
        
        pooledvcov <- pooledvcovout(FullTableFull,betaFE,YYlist,XXlist)
        
        SURbetaFE <- surfeout(FullTableFull,pooledvcov,testvars,nontestvars,
		    YYlist,XXlist)
        
        list2env(wfeandiout(FullTableFull,pooledvcov,testvars,nontestvars,
		    SURbetaFE,YYlist,XXlist), envir = environment())
        
        TrackGroup <- trackgroupout(FullTableFull,pooledvcov,testvars,
		    nontestvars,sigmai,betai,betaWFE,TrackGroup,YYlist,XXlist)
        
        print(counting)
        counting <- counting + 1
        
        Difference <- Initgroups$Group - TrackGroup$Group
        Initgroups <- TrackGroup
        
        if (counting > 15) {
            break
        }
        
    }
	
    Differentstarts <- Initgroups
    FullTableFull$Group <- NULL
    FullTableFull <- merge(FullTableFull, Initgroups, by = c("Vessel"))
	FullTableFull <- FullTableFull[order(FullTableFull$obsnum),]
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    
    betaFE <- indEQFEout(FullTableFull,YYlist,XXlist)
    
    pooledvcov <- pooledvcovout(FullTableFull,betaFE,YYlist,XXlist)
    
    SURbetaFE <- surfeout(FullTableFull,pooledvcov,testvars,nontestvars,YYlist,
	    XXlist)
    
    list2env(wfeandiout(FullTableFull,pooledvcov,testvars,nontestvars,SURbetaFE,
	    YYlist,XXlist), envir = environment())
    
    testDStemp <- testDSout(FullTableFull,pooledvcov,testvars,nontestvars,
	    sigmai,betai,betaWFE,Vessel,YYlist,XXlist)
    
    testDStrack <- testDStemp
    
    return(list(Differentstarts, testDStrack))
    
}
