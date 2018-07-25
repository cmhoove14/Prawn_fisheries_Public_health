# Description  
This folder contains scripts, plots, and simulations for the prawn aquaculture model simulating dynamics of M. rosenbergii and M. volenvenhovenii under repeated stocking amd harvesting events  

## `prawn_aquaculture_mod.R`  
Contains the base model with state variables number of prawns, P, and average prawn length, L; prawn model parameters; sample starting conditions; and a run to equilibrium

## `macrobrachium_aquaculture_data.R`  
Data from Ranjeet and Kurup 2006 on M. rosenbergii stocking trials at different densities used to derive density dependence and marketable biomass from the harvest for the prawn aquaculture model  

## `fit_dens_dep_.R`  
Qualitative fitting of density dependent parameters to data from Ranjeet and Kurup 2006 on M. rosenbergii stocking trials  

## `prawn_aquaculture_sims.R`  
Script to run simulations used to derive eumetric curve (identifies profit-optimal initial stocking density for each species) and then simulate a two year aquaculture run for each species given the profit-optimal conditions identified  

## `prawn_aquaculture_plot.R`  
Takes simulations data and produces manuscript figures 2 and 3 showing the eumetric curve and parameters of the prawn model over an aquqculture cycle

## Other scripts  
Deprecated; mostly old scripts investigating optimization  