#### Full snail-prawn model including epidemiological, predation, and aquaculture components

source('Combined_Model/epi_no_diag_prawn_mod.R')
load("Prawn_aquaculture/aquaculture_sims.RData")

## Set initial values and parameters ###########
years = 10   #number of years to simulate
t.all = c(0:(years*2*365)+366)  #total time vector

#Volenhovenii stocking and harvesting events ###########
nstart.vol = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
               E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
               I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
               Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu, 
               P = opt.vol$P_nought, L = opt.vol$L_nought) #Optimal stocking parameters from volenhovenii ROI optimization

  harvest.t.vol = opt.vol$h.t           #Optimal harvest time from volenhovenii optimization
  n.harvest.vol = floor(((years*365)+366)/harvest.t.vol)
  stocks.vol = data.frame(var = rep(c('P', 'L'), n.harvest.vol),
                          time = rep(366 + harvest.t.vol*c(0:(n.harvest.vol-1)), each = 2),
                          value = rep(c(nstart.vol['P'], nstart.vol['L']), n.harvest.vol),
                          method = rep('rep', n.harvest.vol*2))
    vol.harvests = stocks.vol[seq(1, nrow(stocks.vol), 2), 2]

  par.all.vol = c(par.aqua, par.snails,
                  a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  
#Rosenbergii stocking and harvesting events ###########
nstart.ros = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
               E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
               I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
               Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu, 
               P = opt.ros$P_nought, L = opt.ros$L_nought) #Optimal stocking parameters from volenhovenii ROI optimization

  harvest.t.ros = opt.ros$h.t
  n.harvest.ros = floor(((years*365)+366)/harvest.t.ros)
  stocks.ros = data.frame(var = rep(c('P', 'L'), n.harvest.ros),
                          time = rep(366 + harvest.t.ros*c(0:(n.harvest.ros-1)), each = 2),
                          value = rep(c(nstart.ros['P'], nstart.ros['L']), n.harvest.ros),
                          method = rep('rep', n.harvest.ros*2))
  ros.harvests = stocks.ros[seq(1, nrow(stocks.ros), 2), 2]
  
  par.all.ros = c(par.aqua, par.snails,
                  a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  
  par.all.ros['k'] = 0.0104333333  # alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333
  

#allow for one year run in to display epi model at eqbm ###########
t.yr1 = c(0:365) 
nstart.yr1 = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
               E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
               I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
               Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu)

yr1 = as.data.frame(ode(nstart.yr1,t.yr1,snail_epi_allvh,par.snails))

  yr1$P = 0   #Add prawn variable to integrate with prawn stocking sims
  yr1$L = 0   #Add prawn variable to integrate with prawn stocking sims
  yr1$W = cov*yr1$Wt + (1-cov)*yr1$Wu
  yr1$prev = pnbinom(2, size = 0.2, mu = yr1$W, lower.tail = FALSE)   
  yr1$S.t = (yr1$S1 + yr1$S2 + yr1$S3) / area        # density susceptible snails
  yr1$E.t = (yr1$E1 + yr1$E2 + yr1$E3) / area        # density exposed snails 
  yr1$I.t = (yr1$I2 + yr1$I3 ) / area                # density infected snails
  yr1$N.t = (yr1$S.t + yr1$E.t + yr1$I.t)            # density snails
  yr1$t.1 = (yr1$S1 + yr1$E1) / area                    # density snails of size class 1
  yr1$t.2 = (yr1$S2 + yr1$E2 + yr1$I2) / area      # density snails of size class 2
  yr1$t.3 = (yr1$S3 + yr1$E3 + yr1$I3) / area      # density snails of size class 3

#Simulate annual MDA for n years ########
#mda events
  mdas = data.frame(var = rep('Wt', years),
                    time = 365*c(1:years)+1,
                    value = rep((1-cov*eff), years),
                    method = rep('multiply', years))
  
  sim.mda = as.data.frame(ode(nstart.yr1,t.all,snail_epi_allvh,par.snails,
                          events = list(data = mdas)))
  
  sim.mda$P = 0   #Add prawn variable to integrate with prawn stocking sims
  sim.mda$L = 0   #Add prawn variable to integrate with prawn stocking sims
  sim.mda$W = cov*sim.mda$Wt + (1-cov)*sim.mda$Wu
  sim.mda$prev = pnbinom(2, size = 0.2, mu = sim.mda$W, lower.tail = FALSE)   
  sim.mda$S.t = (sim.mda$S1 + sim.mda$S2 + sim.mda$S3) / area        # density susceptible snails
  sim.mda$E.t = (sim.mda$E1 + sim.mda$E2 + sim.mda$E3) / area        # density exposed snails 
  sim.mda$I.t = (sim.mda$I2 + sim.mda$I3 ) / area                    # density infected snails
  sim.mda$N.t = (sim.mda$S.t + sim.mda$E.t + sim.mda$I.t)            # density snails
  sim.mda$t.1 = (sim.mda$S1 + sim.mda$E1) / area                    # density snails of size class 1
  sim.mda$t.2 = (sim.mda$S2 + sim.mda$E2 + sim.mda$I2) / area      # density snails of size class 2
  sim.mda$t.3 = (sim.mda$S3 + sim.mda$E3 + sim.mda$I3) / area      # density snails of size class 3
  
  sim.mda = rbind(yr1, sim.mda)
    
#Simulate M. volenhovenii stocking for n years ########
sim.vol = as.data.frame(ode(nstart.vol,t.all,snail_prawn_model,par.all.vol,
                              events = list(data = stocks.vol)))
  
  sim.vol$W = cov*sim.vol$Wt + (1-cov)*sim.vol$Wu
  sim.vol$prev = pnbinom(2, size = 0.2, mu = sim.vol$W, lower.tail = FALSE)   
  sim.vol$S.t = (sim.vol$S1 + sim.vol$S2 + sim.vol$S3) / area        # density susceptible snails
  sim.vol$E.t = (sim.vol$E1 + sim.vol$E2 + sim.vol$E3) / area        # density exposed snails 
  sim.vol$I.t = (sim.vol$I2 + sim.vol$I3 ) / area                    # density infected snails
  sim.vol$N.t = (sim.vol$S.t + sim.vol$E.t + sim.vol$I.t)            # density snails
  sim.vol$t.1 = (sim.vol$S1 + sim.vol$E1) / area                          # density snails of size class 1
  sim.vol$t.2 = (sim.vol$S2 + sim.vol$E2 + sim.vol$I2) / area      # density snails of size class 2
  sim.vol$t.3 = (sim.vol$S3 + sim.vol$E3 + sim.vol$I3) / area      # density snails of size class 3
  
  sim.vol = rbind(yr1, sim.vol)
  
#Simulate M. rosenbergii stocking for n years ########
  sim.ros = as.data.frame(ode(nstart.ros,t.all,snail_prawn_model,par.all.ros,
                              events = list(data = stocks.ros)))
  
  sim.ros$W = cov*sim.ros$Wt + (1-cov)*sim.ros$Wu
  sim.ros$prev = pnbinom(2, size = 0.2, mu = sim.ros$W, lower.tail = FALSE)   
  sim.ros$S.t = (sim.ros$S1 + sim.ros$S2 + sim.ros$S3) / area        # density susceptible snails
  sim.ros$E.t = (sim.ros$E1 + sim.ros$E2 + sim.ros$E3) / area        # density exposed snails 
  sim.ros$I.t = (sim.ros$I2 + sim.ros$I3 ) / area                    # density infected snails
  sim.ros$N.t = (sim.ros$S.t + sim.ros$E.t + sim.ros$I.t)            # density snails
  sim.ros$t.1 = (sim.ros$S1 + sim.ros$E1) / area                          # density snails of size class 1
  sim.ros$t.2 = (sim.ros$S2 + sim.ros$E2 + sim.ros$I2) / area      # density snails of size class 2
  sim.ros$t.3 = (sim.ros$S3 + sim.ros$E3 + sim.ros$I3) / area      # density snails of size class 3
  
  sim.ros = rbind(yr1, sim.ros)

#Simulate annual mda along with volenhovenii stocking #############
mda.vol = rbind(mdas, stocks.vol)
  mda.vol = mda.vol[order(mda.vol$time),]
  
  sim.mda.vol = as.data.frame(ode(nstart.vol,t.all,snail_prawn_model,par.all.vol,
                              events = list(data = mda.vol)))
  
  sim.mda.vol$W = cov*sim.mda.vol$Wt + (1-cov)*sim.mda.vol$Wu
  sim.mda.vol$prev = pnbinom(2, size = 0.2, mu = sim.mda.vol$W, lower.tail = FALSE)   
  sim.mda.vol$S.t = (sim.mda.vol$S1 + sim.mda.vol$S2 + sim.mda.vol$S3) / area        # density susceptible snails
  sim.mda.vol$E.t = (sim.mda.vol$E1 + sim.mda.vol$E2 + sim.mda.vol$E3) / area        # density exposed snails 
  sim.mda.vol$I.t = (sim.mda.vol$I2 + sim.mda.vol$I3 ) / area                    # density infected snails
  sim.mda.vol$N.t = (sim.mda.vol$S.t + sim.mda.vol$E.t + sim.mda.vol$I.t)            # density snails
  sim.mda.vol$t.1 = (sim.mda.vol$S1 + sim.mda.vol$E1) / area                          # density snails of size class 1
  sim.mda.vol$t.2 = (sim.mda.vol$S2 + sim.mda.vol$E2 + sim.mda.vol$I2) / area      # density snails of size class 2
  sim.mda.vol$t.3 = (sim.mda.vol$S3 + sim.mda.vol$E3 + sim.mda.vol$I3) / area      # density snails of size class 3
  
  sim.mda.vol = rbind(yr1, sim.mda.vol)

#Simulate annual mda along with rosenbergii stocking #############
mda.ros = rbind(mdas, stocks.ros)
  mda.ros = mda.ros[order(mda.ros$time),]
  
  sim.mda.ros = as.data.frame(ode(nstart.ros,t.all,snail_prawn_model,par.all.ros,
                                  events = list(data = mda.ros)))
  
  sim.mda.ros$W = cov*sim.mda.ros$Wt + (1-cov)*sim.mda.ros$Wu
  sim.mda.ros$prev = pnbinom(2, size = 0.2, mu = sim.mda.ros$W, lower.tail = FALSE)   
  sim.mda.ros$S.t = (sim.mda.ros$S1 + sim.mda.ros$S2 + sim.mda.ros$S3) / area        # density susceptible snails
  sim.mda.ros$E.t = (sim.mda.ros$E1 + sim.mda.ros$E2 + sim.mda.ros$E3) / area        # density exposed snails 
  sim.mda.ros$I.t = (sim.mda.ros$I2 + sim.mda.ros$I3 ) / area                    # density infected snails
  sim.mda.ros$N.t = (sim.mda.ros$S.t + sim.mda.ros$E.t + sim.mda.ros$I.t)            # density snails
  sim.mda.ros$t.1 = (sim.mda.ros$S1 + sim.mda.ros$E1) / area                          # density snails of size class 1
  sim.mda.ros$t.2 = (sim.mda.ros$S2 + sim.mda.ros$E2 + sim.mda.ros$I2) / area      # density snails of size class 2
  sim.mda.ros$t.3 = (sim.mda.ros$S3 + sim.mda.ros$E3 + sim.mda.ros$I3) / area      # density snails of size class 3
  
  sim.mda.ros = rbind(yr1, sim.mda.ros)
    
  
save.image(file = "Combined_Model/epi_no_diag_prawn_sims_n-2.RData")  