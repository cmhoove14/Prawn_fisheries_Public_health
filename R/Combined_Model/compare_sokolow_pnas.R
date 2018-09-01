#Script to check elimination thresholds are consistent with those presented in Sokolow et al PNAS 2015

source("Data/aquaculture_parameters.R")
source("Data/snail_epi_parameters.R")
source("Data/combined_parameters.R")

source("R/Combined_Model/epi_no_diag_prawn_immigration_mod.R")
source("R/Epi_Model/snail_epi_mod_no_diag_immigration.R")

require(deSolve)
require(tidyverse)

# First run epi model to eqbm by itself
epi_eqbm <- as.data.frame(ode(nstart.sn,
                              seq(1,365*30,30),
                              snail_epi_allvh_imm,
                              par.snails.imm))

nstart.eqbm <- epi_eqbm %>% filter(time == max(time)) %>% select(S1:Wu) %>% unlist()

#First turn off mortality and growth of prawns 
par.all = c(par.aqua, par.snails.imm, par.epi_prawn)
  par.all["muP"] <- 0
  par.all["om"] <- 0
  par.all["k"] <- 0
  par.all["nfr"] <- 1   #Holling's type II as in Sokolow et al
  par.all["eps"] <- 10  #Baseline, but make sure to assign here for reproducibility as it may change based on the results of this script
  
# Function to run model to equilibrium with constant prawn density 
sim_constant_prawn <- function(P_dens, var){
  nstart <- nstart.eqbm
  nstart["P"] <- P_dens * area
  nstart["L"] <- 150   #Adult prawns
  
  run_eqbm <- as.data.frame(ode(nstart,
                                seq(1,365*30,30),
                                snail_prawn_model_imm,
                                par.all)) %>% 
    mutate(W = cvrg*Wt + (1-cvrg)*Wu,
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t)) # density snails of size class 3
  
  return(run_eqbm %>% filter(time == max(time)) %>% pull(!!enquo(var)))
}

constant_prawn_sims <- data.frame(P_dens = seq(0,1,0.025)) %>% 
  mutate(W = map_dbl(P_dens, sim_constant_prawn, W),
         N_dens = map_dbl(P_dens, sim_constant_prawn, N.t))

constant_prawn_sims %>% ggplot(aes(x = P_dens, y = N_dens)) +
  geom_line(size = 1.2) + theme_bw()

#Prawns are packing way too much punch, increase the wildland consumption rate penalty
  par.all["eps"] <- 50

constant_prawn_sims2 <- data.frame(P_dens = seq(0,1,0.025)) %>% 
  mutate(W = map_dbl(P_dens, sim_constant_prawn, W),
         N_dens = map_dbl(P_dens, sim_constant_prawn, N.t))

constant_prawn_sims2 %>% ggplot(aes(x = P_dens, y = N_dens)) +
  geom_line(size = 1.2) + theme_bw()

constant_prawn_sims2 %>% ggplot(aes(x = P_dens, y = W)) +
  geom_line(size = 1.2, col = "purple") + theme_bw()

