# Description  
Schisto epidemiological model with infection and size stratified snail population that gives rise to adult worms harbored in definitive human hosts. Simulations and plots to explore dynamics

## `snail_epi_mod_no_diag_immigration`  
Epidemiological model simulating infection and growth dynamics across three infection states and three size classes with all vertical and horizontal transitions possible, but no diagonal transitions based on the assumption that it is unlikely on a daily time scale for a particular snail to transition in both growth and infection status. Also includes functionality for migration with emigration away from the site equal to immigration to the site  

## `snail_epi_life_expectancy`  
Simulations and tweaks to snail mortality parameters to get life expectancy and size/infection dynamics right  

## `snail_spi_sim`  
Function to run model to euqilibrium given parameter set and starting conditions and return worm burden and density of each snail infection class  

## `Archive`  
Older versions of the model and other scripts for simulation and plotting
