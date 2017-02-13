<!-- README.md is generated from README.Rmd. Please edit that file -->
BDAepimodel
===========

Efficient Data Augmentation for Fitting Stochastic Epidemic Models to Prevalence Data
-------------------------------------------------------------------------------------

The **BDAepimodel** R package implements a Bayesian data augmentation algorithm for fitting stochastic epidemic models with arbitrary dynamics to disease prevalence data. The implementation is reasonably fast for progressive, non-recurrent stochastic epidemic models in populations as large as several thousand individuals, but is not yet optimized for large population inference. This repository and package exist to expose and archive the code used in the 2016 paper, "Efficient Data Augmentation for Fitting Stochastic Epidemic Models to Prevalence Data" by Fintzi, Wakefield, and Minin. Future work and extensions will be implemented in the 'stemr' package (<https://github.com/fintzij/stemr>).

Getting started
---------------

This package may be installed directly from GitHub using the **devtools** package:

    library(devtools)
    install_github("fintzij/ECctmc",build_vignettes=TRUE) # required for simulating sample paths for endpoint conditioned CTMCs
    install_github("fintzij/BDAepimodel",build_vignettes=TRUE) 
    library(BDAepimodel)

There are several vignettes included in this package. The first diagrams the basic functionality of the package and is intended to serve as a how-to for setting up MCMC using our code. The others contain walk-throughs of the simulations presented in the original paper, and of the analysis of the influenza outbreak data that is presented in the paper. These are accessible as follows:

    vignette("BDAepimodel") 
    vignette("SIR_SEIR_SIRS")
    vignette("model_misspecification")
    vignette("popsize_misspecification")
    vignette("prior_effect")
    vignette("bbs_influenza")
    
We also provide vignettes that demonstrate how to fit the models in the first and second simulations, along with the boarding school examples, using the pomp package. These are accessible as follows:

    vignette("SIR_pomp_simulation")
    vignette("SEIR_pomp_simulation")
    vignette("SIRS_pomp_simulation")
    vignette("model_misspecification_pomp")
    vignette("bbs_pomp_binom_vignette")
    vignette("bbs_pomp_negbinom_vignette")
