#Function to simulate aquaculture cycle with given parameters and starting conditions and return optimal management conditions 
sim_aqua_eum <- function(a.p, b.p, gam, muP, d, om, k, linf, k.ros, c, fc, p, delta, L_nought, P_nought, species){
  
  #Set starting parameters based on stocking density  
    sim_start = c(P = P_nought, L = L_nought)
  
  #Get marketable fraction of harvest as function of stocking density    
    market_frac = predict(eta.lm, newdata = data.frame(dens = P_nought/area))
    
  # Set parameters from row of eum_sim
    # Set growth rate based on species  
      if(species == "M. vollenhovenii"){
        par.aqua['k'] = k
      } else {
        par.aqua['k'] = k.ros
      }

      par.aqua['a.p'] <- a.p
      par.aqua['b.p'] <- b.p
      par.aqua['gam'] <- gam
      par.aqua['muP'] <- muP
      par.aqua['d'] <- d
      par.aqua['om'] <- om
      par.aqua['linf'] <- linf
      
      cost <- c
      fixed_cost <- fc
      price <- p
      delta <- delta
    
  #Simulate for two years  
    sim_df <- as.data.frame(ode(sim_start,t.p,prawn_biomass,par.aqua)) %>% 
      mutate(p_mass = 10^(a.p+b.p*log10(L/10)),     
             total_mass = (P*p_mass) /1000,    
             harvest_mass = total_mass * market_frac,
             revenue = harvest_mass*price,
             profit = revenue*exp(-delta*time) - cost*P_nought - fixed_cost,
             roi = profit / (cost*P_nought),
             n_harvest = floor((years*365)/time),
             cum_profits = n_harvest*profit,
             Species = species) %>% 
      filter(p_mass >= 30)
    
    op_mgmt <- sim_df %>% filter(cum_profits == max(cum_profits))
    
  return(op_mgmt) 

}

#Function to simulate aquaculture cycle through time 
sim_aqua_time <- function(P_nought, species = "M. vollenhovenii"){
  
#Set starting parameters based on stocking density  
  sim_start = c(P = P_nought, L = 40)

#Get marketable fraction of harvest as function of stocking density    
  market_frac = predict(eta.lm, newdata = data.frame(dens = P_nought/area))
  
#Set growth rate based on species  
  if(species == "M. vollenhovenii"){
    par.aqua['k'] = k.vol
  } else {
    par.aqua['k'] = k.ros
  }
  
#Simulate for two years  
  sim_df <- as.data.frame(ode(sim_start,t.p,prawn_biomass,par.aqua)) %>% 
      mutate(P_dens = P/area,
             P_start = P_nought / area,
             p_mass = 10^(par.aqua["a.p"]+par.aqua["b.p"]*log10(L/10)),     
             total_mass = (P*p_mass) /1000,    
             harvest_mass = total_mass * market_frac,
             revenue = harvest_mass*price,
             profit = (revenue*exp(-delta*time)) - (cost*P_nought) - fixed_cost,
             roi = profit / (cost*P_nought),
             n_harvest = floor((years*365)/time),
             cum_profits = n_harvest*profit,
             Species = species)

  return(as.matrix(sim_df))  

}

