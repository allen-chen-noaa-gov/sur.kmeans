kmeans_wrapper <- function(FullTableFull,YYlist,XXlist,Vessel,testvars,
    nontestvars,nseeds,startseeds,ntreat=1,runparallel=FALSE,glimit=3) {
    #' Outside wrapper for full K-means algorithm
    #'
    #' Outside wrapper
    #'
    #' @param FullTableFull data
    #' @return Converged groupings and associated statistical tests
    #' @export
    #' @examples
    #'
	
	testcheck = 0
    FullTableFull$obsnum <- 1:dim(FullTableFull)[1]

	while (length(Vessel)/ntreat >= glimit & testcheck == 0) {
	
	    if(ntreat==1) {
		
            betaFE <- indEQFEout(FullTableFull,YYlist,XXlist)

            pooledvcov <- pooledvcovout(FullTableFull,betaFE,YYlist,XXlist)

            SURbetaFE <- surfeout(FullTableFull,pooledvcov,testvars,nontestvars,
		        YYlist,XXlist)

            #wfe, individual beta, and individual sigma (for dispersion test)
            list2env(wfeandiout(FullTableFull,pooledvcov,testvars,nontestvars,
		        SURbetaFE,YYlist,XXlist),envir=environment())

            #dispersion test
            testDS <- testDSout(FullTableFull,pooledvcov,testvars,nontestvars,
			    sigmai,betai,betaWFE,Vessel,YYlist,XXlist)

            print(testDS)
			
			#two-sided standard normal
            if (2 * pnorm(-abs(testDS[[1]])) < 0.05) {
			
                print("REJECT NULL")   
				testcheck <- 0
				
            } else {
			
                print("FAIL TO REJECT NULL")
				testcheck <- 1

            }
		
		    print(Sys.time())

            ntreat <- ntreat+1

        } else {
	        
		    res <- kmeansout(FullTableFull,ntreat,nseeds,Vessel,startseeds,
			    testvars,nontestvars,runparallel,YYlist,XXlist)

            list2env(res,envir=environment())

            ntreat <- ntreat+1
        
		    print(Sys.time())
		
        }

	}
	
	if (length(Vessel)/ntreat < glimit) {
        print("Too few individuals per group")
    } else if (testcheck == 1) {
	    print("Converged")
	}
	
	return(list(index = index, FinalGrouping = FinalGrouping,
	    testDStrack = testDStrack, Differentstarts = Differentstarts))
	
}
