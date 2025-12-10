# Bathymetry

This project operates as a consistent interface for [aegis.* projects](https://github.com/jae0/) to bathymetry data and modelled predictions. 

Currently, it is used for bathymetric data and modelled predictions for the North-West Atlantic near the focal areal of the aegis.* projects, but this can be readily altered for your area of interest. You can use the examples as a template or use your own pre-existing data to format to the expectations of the interpolation methods used.

Here the primary methods used for spatial interpolation include:

1. [stmv - Space-time models of variability](https://github.com/jae0/stmv) - [a flexible framework for spatial (and spatiotemporal) interpolation](https://github.com/jae0/stmv/blob/master/docs/stmvMethods.pdf) using continuous autocorrelation models that are discretized to an arbitary resolution (raster). It acts as a front-end to interpolation engines such as Gaussian random fields (via INLA, RandomFields), Kernal-based methods (via FFTW), and classical kriging with external drift (gstat) that has a moving window of operation to account for *non-stationary* spatial processes. It uses mgcv or GLM to model the "external drift" covariate effects and the spatial (and spatiotemporal) models to account for the random autocorrelation processes. It is, therefore, an ad-hoc method that however [works surprisingly well](https://github.com/jae0/aegis.bathymetry/blob/master/inst/scripts/02_bathymetry_stmv.R), albeit, slowly. However, as computational power increases, Expectation Maximization or Gibbs sampling would provide a more rigourous foundation. A faster version on a [coarser resolution can be found here](https://github.com/jae0/aegis.bathymetry/blob/master/inst/scripts/99_bathymetry_stmv_example.R).  

2. [carstm - Conditional autoregressive space-time models](https://github.com/jae0/carstm) - an implementation of CAR (Conditional AutoRegressive) and BYM (Besag-York_Mollie) approach to areal unit models. It is a convenience front end to INLA that inter-operates with aegis.* projects in a coherent manner to permit more expressive models. [Examples of how to use it are shown here](https://github.com/jae0/carstm/blob/master/inst/scripts/example_temperature_carstm.R) where ocean bottom temperatures are modelled across space (and time). The method for bathymetry itself is [found here](https://github.com/jae0/aegis.bathymetry/blob/master/inst/scripts/03_bathymetry_carstm.R).




