SUR Kmeans
=========
---

If you run into problems you can contact allen.chen@noaa.gov. Latent class inference with k-means and homogeneity testing implemented for a Seemingly Unrelated Regressions
framework: uses Lin and Ng's (2012 Journal of Econometric Methods) k-means algorithm in conjunction with Pesaran and Yamagata's (2008 Journal of Econometrics) parameter
homoegeneity test to uncover latent classes, implemented here for a Seemingly Unrelated Regressions framework.

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
