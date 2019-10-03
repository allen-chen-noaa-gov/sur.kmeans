kmeansout <- function(FullTableFull,ntreat,nseeds,Vessel,startseeds,testvars,
    nontestvars,runparallel = FALSE,YYlist,XXlist) {
    #' K-means algorithm to find latent groupings
    #'
    #' Keeps regrouping until convergence
    #'
    #' @param FullTableFull data
    #' @return Converged groupings and associated statistical tests
    #' @export
    #' @examples
    #'
    
	if(runparallel==TRUE) {
	    cores <- tryCatch({
		    if(!requireNamespace("foreach", quietly = TRUE)) {
                stop("The foreach package is required for parallel scenario
				    simulations.", call. = FALSE)
            }
			cores <- foreach::getDoParWorkers()
			if(cores == 1) {
                warning("You have only registered one core. Reverting to
				    non-parallel processing.")
            }
			cores
        }, error = function(e) {
            1			
        })
        if(cores == 1) {
		    runparallel <- FALSE
	    }
		# Code from ss3sim setup_parallel function
    }
	
	if(runparallel==TRUE) {
        par_res <- foreach(nn = startseeds:nseeds,
	        .packages = c("sur.kmeans")) %dopar% kmeans_inner(nn,FullTableFull,
		    ntreat,Vessel,testvars,nontestvars,YYlist,XXlist)
	} else {
        par_res <- lapply(startseeds:nseeds,kmeans_inner,
		    FullTableFull=FullTableFull,ntreat=ntreat,Vessel=Vessel,
			testvars=testvars,nontestvars=nontestvars,YYlist,XXlist)
    }
    
    Differentstarts <- c(rep(list(NA), startseeds - 1),
	    lapply(par_res, "[[", 1))
    testDStrack <- c(rep(list(Inf), startseeds - 1), lapply(par_res, "[[", 2))
    
    meantests <- list()
    for (j in 1:nseeds) {
        meantests[[j]] <- sum(abs(as.numeric(testDStrack[[j]])))
    }
    # start from 1 not `startseeds` because we've replaced missing values with
	    # NAs and Infs
    
    test <- Reduce("min", meantests)
    index <- match(test, meantests)
    
    FinalGrouping <- Differentstarts[[index]]
    
    FullTableFull$Group <- NULL
    FullTableFull <- merge(FullTableFull, FinalGrouping, by = c("Vessel"))
	FullTableFull <- FullTableFull[order(FullTableFull$obsnum),]
    FullTableFull$Group <- as.factor(FullTableFull$Group)
    
    testDS <- testDStrack[[index]]
    
    print(testDS)
	
    #two-sided standard normal
    if (any(2 * pnorm(-abs(unlist(testDS))) < 0.05)) {
	
        print("REJECT NULL")
        testcheck <- 0
	
    } else {
	
        print("FAIL TO REJECT NULL")
        testcheck <- 1
	
    }
    
    return(list(testcheck = testcheck, index = index,
	    FinalGrouping = FinalGrouping, testDStrack = testDStrack,
		Differentstarts = Differentstarts))
	
}
