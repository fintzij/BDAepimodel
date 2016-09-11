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
    install_github("fintzij/BDAepimodel",build_vignettes=TRUE) 
    library(BDAepimodel)
    vignette("BDAepimodel")
