# Function to simulate model for prawn only and combined prawn and mda interventions ################
# Returns data frame of full simulation through time, used for plots comparing simulations through time
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
           W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu,
           prev = get_prev(parameters["phi"], W),
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
                                 events = list(data = stock.events),
                                 atol = 1e-6, rtol = 1e-6, method = "radau")) %>% 
    mutate(W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu,
           prev = get_prev(parameters["phi"], W),
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
           W = parameters["cvrg"]*Wt + (1-parameters["cvrg"])*Wu,
           prev = get_prev(parameters["phi"], W),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3
  
#Data frame of full simulation through runup time, intervene time, and post intervention time  
  sim.fin <- rbind(runup[-max(t.runup),], intervene[-which(intervene$time == max(t.intervene)),], post)
  
  return(sim.fin)
}

# Function to simulate prawn intervention model with a particular parameter set and return outcomes of interest ###########
# Used in global sensitivity analysis
prawn_intervention_sim <- function(pars){
  
  opt_df = expand.grid.df(data.frame(pars[c("a.p", "b.p", "gam", "muP", "d", "om", 
                                            "k", "linf", "k.ros",
                                            "c", "fc", "p", "eta", "delta")]), 
                          data.frame(L_nought = 40, P_nought = seq(500, 7500, 100))) 

#Get mgmt strategy that maximizes cumulative ten year profits
opt_stock <- opt_df %>% 
  bind_cols(pmap_df(., sim_aqua_eum, species = "M. rosenbergii")) %>% 
  filter(cum_profits == max(cum_profits))

# Generate data fram of stocking/harvesting events
  harvest_time <- opt_stock$time
  n.harvest = opt_stock$n_harvest
  stocks = data.frame(var = rep(c('P', 'L'), n.harvest),
                          time = rep(harvest_time*c(0:(n.harvest-1))+1, each = 2),
                          value = rep(c(opt_stock$P, opt_df$L_nought[1]), n.harvest),
                          method = rep('rep', (n.harvest)*2))
  
# Run epi model to equilibrium 
  epi_sim <- as.data.frame(ode(nstart,
                               seq(1,365*100,50),
                               snail_epi_allvh_imm,
                               pars))
  
  epi_eqbm <- epi_sim %>% filter(time == max(time)) 
    full_eqbm <- epi_eqbm %>% select(-time) %>% unlist() 
    
#Add in prawn info and migration info to starting conditions/parameter set    
  full_eqbm["P"] <- opt_stock$P
  full_eqbm["L"] <- opt_df$L_nought[1] 
    
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
                                        events = list(data = stocks),
                                        atol = 1e-6, rtol = 1e-6, method = "radau")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["phi"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu,
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area)
  
  return(intervention_sim %>% 
           select(W, DALYs_W, N.t, I.t) %>% 
           rbind(c(sum(.$W), sum(.$DALYs_W), sum(.$N.t), sum(.$I.t))) %>% 
           slice(max(c(1:(years*365))+1)))
    
}

# Function to compare different interventions across parameter sets ##############
# used for boxplot comparing cumulative impact of each intervention 
compare_interventions <- function(pars, years){
  
  opt_df = expand.grid.df(data.frame(t(pars[c("a.p", "b.p", "gam", "muP", "d", "om", 
                                            "k", "linf", "k.ros",
                                            "c", "fc", "p", "delta")])), 
                          data.frame(L_nought = 40, P_nought = seq(500, 7500, 100))) 

#Get mgmt strategy that maximizes cumulative ten year profits
opt_stock <- opt_df %>% 
  bind_cols(pmap_df(., sim_aqua_eum, species = "M. rosenbergii")) %>% 
  filter(cum_profits == max(cum_profits))

# Generate data fram of stocking/harvesting events
  harvest_time <- opt_stock$time
  n.harvest = opt_stock$n_harvest
  stocks = data.frame(var = rep(c('P', 'L'), n.harvest),
                          time = rep(harvest_time*c(0:(n.harvest-1))+1, each = 2),
                          value = rep(c(opt_stock$P, opt_df$L_nought[1]), n.harvest),
                          method = rep('rep', (n.harvest)*2))

#mda events #######
  mdas = data.frame(var = rep('Wt', years),
                    time = 365*c(0:(years-1))+1,
                    value = rep((1-pars["eff"]), years),
                    method = rep('multiply', years))
  
#Combined prawn stocking and mda events ################
  mda.prawn = rbind(mdas, stocks)
  mda.prawn = mda.prawn[order(mda.prawn$time),]

# Run epi model to equilibrium  ##############
  epi_sim <- as.data.frame(ode(nstart,
                               seq(1,365*100,50),
                               snail_epi_allvh_imm,
                               pars))
  
  epi_eqbm <- epi_sim %>% filter(time == max(time)) 
    full_eqbm <- epi_eqbm %>% select(-time) %>% unlist() 
  
      
#Add in prawn info and migration info to starting conditions/parameter set    
  full_eqbm["P"] <- opt_stock$P
  full_eqbm["L"] <- opt_df$L_nought[1] 
    
  pars["siteS1"] <- epi_eqbm %>% pull(S1)
  pars["siteE1"] <- epi_eqbm %>% pull(E1)
  pars["siteI1"] <- epi_eqbm %>% pull(I1)
  pars["siteS2"] <- epi_eqbm %>% pull(S2)
  pars["siteE2"] <- epi_eqbm %>% pull(E2)
  pars["siteI2"] <- epi_eqbm %>% pull(I2)
  pars["siteS3"] <- epi_eqbm %>% pull(S3)
  pars["siteE3"] <- epi_eqbm %>% pull(E3)
  pars["siteI3"] <- epi_eqbm %>% pull(I3)
  
  
#Simulate no intervention   ##########
    sim.nothing = as.data.frame(ode(epi_eqbm %>% select(-time) %>% unlist(),
                                    c(1:(years*365)),
                                    snail_epi_allvh_imm,
                                    pars)) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["phi"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)#,
           #S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           #E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           #I.t = (I1 + I2 + I3 ) / area, # density infected snails
           #N.t = (S.t + E.t + I.t), # total density snails
           #t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           #t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           #t.3 = (S3 + E3 + I3) / area)
  
    nothing_DALYs <- sum(sim.nothing %>% pull(DALYs_W))

#Simulate annual MDA for n years ########
  sim.mda = as.data.frame(ode(epi_eqbm %>% select(-time) %>% unlist(),
                              c(1:(years*365)),
                              snail_epi_allvh_imm,
                              pars,
                              events = list(data = mdas),
                              atol = 1e-6, rtol = 1e-6, method = "radau")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["phi"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)#,
           #S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           #E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           #I.t = (I1 + I2 + I3 ) / area, # density infected snails
           #N.t = (S.t + E.t + I.t), # total density snails
           #t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           #t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           #t.3 = (S3 + E3 + I3) / area)
    
    mda_DALYs <- sum(sim.mda %>% pull(DALYs_W))

#Simulate prawn intervention under optimal management for n years    
  prawn_sim <- as.data.frame(ode(full_eqbm,
                                 c(1:(years*365)),
                                 snail_prawn_model_imm,
                                 pars,
                                 events = list(data = stocks),
                                 atol = 1e-6, rtol = 1e-6, method = "radau")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["phi"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)#,
           #S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           #E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           #I.t = (I1 + I2 + I3 ) / area, # density infected snails
           #N.t = (S.t + E.t + I.t), # total density snails
           #t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           #t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           #t.3 = (S3 + E3 + I3) / area)
  
  prawn_DALYs <- sum(prawn_sim %>% pull(DALYs_W))

#Simulate combined intervention under optimal management and annual MDA for n years    
  cmbnd_sim <- as.data.frame(ode(full_eqbm,
                                 c(1:(years*365)),
                                 snail_prawn_model_imm,
                                 pars,
                                 events = list(data = mda.prawn),
                                 atol = 1e-6, rtol = 1e-6, method = "radau")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           prev = get_prev(pars["phi"], W),
           DALYs_Wt = map_dbl(Wt, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * pars["cvrg"])),
           DALYs_Wu = map_dbl(Wu, est_dalys, 
                              pars["phi"], pars["weight_lo"], pars["weight_hi"],
                              pars["epmL"], round(pars["H"] * (1-pars["cvrg"]))),
           DALYs_W = DALYs_Wt + DALYs_Wu)#,
           #S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           #E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           #I.t = (I1 + I2 + I3 ) / area, # density infected snails
           #N.t = (S.t + E.t + I.t), # total density snails
           #t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           #t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           #t.3 = (S3 + E3 + I3) / area)
  
  cmbnd_DALYs <- sum(cmbnd_sim %>% pull(DALYs_W))
  
  return(data.frame("None" = nothing_DALYs, 
                    "MDA" = mda_DALYs, 
                    "Prawns" = prawn_DALYs, 
                    "Prawns & MDA" = cmbnd_DALYs))
  
}
