SUR Kmeans
=========
---

If you run into problems you can contact allen.chen@noaa.gov

# Installation #
---

To install with vignette:

    > install.packages("devtools")
	> library(devtools)
	
Set the directory you've downloaded the package into:

    > setwd("J:/Fishperson/Directory_containing_sur.kmeans")

No vignettes to install with yet:

    > install("sur.kmeans", build_vignettes = FALSE)
	
Check out documentation:

    > help(package="sur.kmeans")
	
Need to add simulated data for test usage, then run wrapper function with 
    variable names and data:

    > kmeans_wrapper(FullTable, YYlist, XXlist, Vessel, testvars, nontestvars,
	    nseeds, startseeds, runparallel=TRUE)
