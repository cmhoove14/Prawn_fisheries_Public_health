source("../Data/Ranjeet_Kurup_06_data.R")

  eta.lm = lm(marketable ~ dens, data = rk06)


#Small function to estimate cumulative discounted profits given number of cycles, profit per cycle, discount rate, and time per cycle
get_cum_profits <- function(n, Pi, delta, Time){
  if(is.na(Pi)){
    
    return(NA)
    
  } else {
    
    return(sum(sapply(c(1:n), function(n) Pi*exp(-delta*n*Time))))
    
  }
}

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
      #If average prawn mass is less than 30 g (marketable size) don't estimate monetary outcomes       
             revenue = ifelse(p_mass <30, NA, harvest_mass*price),
             profit = ifelse(p_mass <30, NA, revenue - cost*P_nought - fixed_cost),
             roi = ifelse(p_mass <30, NA, profit / (cost*P_nought)),
             n_harvest = ifelse(p_mass <30, NA, floor((years*365)/time)),
             cum_profits = pmap_dbl(list(n = n_harvest, 
                                         Pi = profit, 
                                         delta = delta, 
                                         Time = time), get_cum_profits),
             Species = species)
    
    op_mgmt <- sim_df %>% filter(cum_profits == max(cum_profits, na.rm = T))
    
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
      #If average prawn mass is less than 30 g (marketable size) don't estimate monetary outcomes       
             revenue = ifelse(p_mass <30, NA, harvest_mass*price),
             profit = ifelse(p_mass <30, NA, revenue - cost*P_nought - fixed_cost),
             roi = ifelse(p_mass <30, NA, profit / (cost*P_nought)),
             n_harvest = ifelse(p_mass <30, NA, floor((years*365)/time)),
             cum_profits = pmap_dbl(list(n = n_harvest, 
                                         Pi = profit, 
                                         delta = delta, 
                                         Time = time), get_cum_profits),
             Species = species)

  return(sim_df)  

}

