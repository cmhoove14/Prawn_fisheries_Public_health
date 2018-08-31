# Description  
This folder contains scripts, plots, and simulations for the prawn aquaculture model simulating dynamics of M. rosenbergii and M. volenvenhovenii under repeated stocking amd harvesting events  

## `prawn_aquaculture_mod.R`  
Contains the base model with state variables number of prawns, P, and average prawn length, L

## `macrobrachium_aquaculture_investigate.R`  
Investigation of data from Ranjeet and Kurup 2006 on M. rosenbergii stocking trials at different densities used to derive density dependence and marketable biomass from the harvest for the prawn aquaculture model as well as different parameterizations of allometric length-weight conversion from Lalrinsanga et al 2012 and Kuris et al 1987  

## `fit_dens_dep_.R`  
Qualitative fitting of density dependent parameters to data from Ranjeet and Kurup 2006 on M. rosenbergii stocking trials  

## `prawn_aquaculture_sim.R`  
Functions to simulate the aquaculture model given different parameters and starting conditions and return either optimal conditions or full simulations across time  

## `Archive`  
Deprecated; old scripts  