prawn_sim <- function(t.runup, t.intervene, t.post,
                      model.runup = snail_epi_allvh_imm,
                      model.intervene = snail_prawn_model_imm,
                      model.post = snail_epi_allvh_imm,
                      nstart, parameters, stock.events){
  
# Simulate model for runup time ###########   
  runup <- as.data.frame(ode(nstart,
                             t.runup,
                             model.runup,
                             parameters)) %>% 
    mutate(P = 0,
           L = 0,
           W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = par.snails.imm['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3  
  
  nstart.intervene <- runup %>% filter(time == max(t.runup)) %>% select(S1:L) %>% unlist()
    nstart.intervene["P"] <- stocks.ros %>% filter(var == "P" & time == 366) %>% pull(value)
    nstart.intervene["L"] <- stocks.ros %>% filter(var == "L" & time == 366) %>% pull(value)
 
# Simulate model for intervention time ###########        
  intervene <- as.data.frame(ode(nstart.intervene,
                                 t.intervene,
                                 model.intervene,
                                 parameters,
                                 events = list(data = stock.events))) %>% 
    mutate(W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = par.snails.imm['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3
    
  nstart.post <- intervene %>% filter(time == max(t.intervene)) %>% select(S1:Wu) %>% unlist()

# Simulate model post intervention ######
  post <- as.data.frame(ode(nstart.post,
                            t.post,
                            model.post,
                            parameters)) %>% 
    mutate(P = 0,
           L = 0,
           W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = par.snails.imm['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3
  
  sim.fin <- rbind(runup[-max(t.runup),], intervene[-which(intervene$time == max(t.intervene)),], post)
  
  return(sim.fin)
}

prawn_intervention_sim <- function(pars){
  
#Run prawn model through two years to identify optimal harvest time
   sim_prawns <- as.data.frame(ode(prawn_start,
                                  c(1:(365*2)),
                                  prawn_biomass,
                                  pars)) %>% 
      mutate(p_mass = 10^(pars["a.p"]+pars["b.p"]*log10(L/10)),     
             total_mass = (P*p_mass) /1000,    
             harvest_mass = total_mass * as.numeric(pars["eta"]),
             revenue = harvest_mass*price,
             profit = revenue*exp(-pars["delta"]*time) - pars["c"]*prawn_start["P"] - pars["fc"],
             n_harvest = floor((years*365)/time),
             cum_profits = n_harvest*profit) %>% 
      filter(p_mass >= 30)
    
    opt_harvest <- sim_prawns %>% filter(cum_profits == max(cum_profits))
    
# Generate data fram of stocking/harvesting events
  harvest_time <- opt_harvest$time
  n.harvest = opt_harvest$n_harvest
  stocks = data.frame(var = rep(c('P', 'L'), n.harvest),
                          time = rep(harvest_time*c(1:n.harvest), each = 2),
                          value = rep(c(opt.ros$P_nought, opt.ros$L_nought), n.harvest),
                          method = rep('rep', (n.harvest)*2))
  
    harvests = stocks[seq(1, nrow(stocks), 2), 2]  

# Run epi model to equilibrium 
  epi_sim <- as.data.frame(ode(nstart,
                               seq(1,365*20,40),
                               snail_epi_allvh_imm,
                               pars))
  
  epi_eqbm <- epi_sim %>% filter(time == max(time)) 
    full_eqbm <- epi_eqbm %>% select(-time) %>% unlist() 
    full_eqbm["P"] <- opt.ros$P_nought 
    full_eqbm["L"] <- opt.ros$L_nought
    
  pars["siteS1"] <- epi_eqbm %>% pull(S1)
  pars["siteE1"] <- epi_eqbm %>% pull(E1)
  pars["siteI1"] <- epi_eqbm %>% pull(I1)
  pars["siteS2"] <- epi_eqbm %>% pull(S2)
  pars["siteE2"] <- epi_eqbm %>% pull(E2)
  pars["siteI2"] <- epi_eqbm %>% pull(I2)
  pars["siteS3"] <- epi_eqbm %>% pull(S3)
  pars["siteE3"] <- epi_eqbm %>% pull(E3)
  pars["siteI3"] <- epi_eqbm %>% pull(I3)
  
# Run model with prawn intervention 
  intervention_sim <- as.data.frame(ode(full_eqbm,
                                        c(1:(years*365)),
                                        snail_prawn_model_imm,
                                        pars,
                                        events = list(data = stocks))) %>% 
    mutate(W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = pars['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) %>% # density snails of size class 3
    filter(time >= ((365*years) - 365)) %>% summarise_all(mean)
  
  return(intervention_sim %>% select(W, N.t, I.t))
    
}
