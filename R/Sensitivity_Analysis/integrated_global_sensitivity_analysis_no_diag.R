#Assess sensitivity of each outcome ############
source('Combined_Model/epi_no_diag_prawn_mod.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
source('Epi_Model/snail_epi_mod_no_diag.R')
source('Sensitivity_Analysis/sensitivity_prep.R')
load('Prawn_aquaculture/aquaculture_sims.RData')

library(sensitivity)
library(parallel)

  
#first get outcomes from prawn aquaculture model ###########
  lhcpars.aqua = lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]
  
  lhcpars.aqua.fill = data.frame(h.t = 0,
                                 h.bm = 0,
                                 h.mbm = 0,
                                 revenue = 0,
                                 profit = 0,
                                 roi = 0)
  
#Run model through with parameter sets from latin hypercube
  for(i in 1:nsims){
    parun<-lhcpars.aqua[i,1:length(names(par.aqua))]
    c = lhcpars.aqua[i,(length(names(par.aqua))+1)]
    p = lhcpars.aqua[i,(length(names(par.aqua))+2)]
    eta = lhcpars.aqua[i,(length(names(par.aqua))+3)]
      
    output = as.data.frame(ode(nstart.r, t.p, prawn_biomass, parun))
    output$B = (par.aqua['a.p']/10)*(output$L/10)^par.aqua['b.p']       
    output$Bt = output$B*output$P      
    
    lhcpars.aqua.fill[i,1] = output$time[output$Bt==max(output$Bt)]   #harvest time
    lhcpars.aqua.fill[i,2] = max(output$Bt)                       #harvest total biomass
    lhcpars.aqua.fill[i,3] = max(output$Bt) * eta                 #harvest marketable biomass
    lhcpars.aqua.fill[i,4] = p*(max(output$Bt)/1000)* eta               # Raw Profit (revenue)
    lhcpars.aqua.fill[i,5] = p*(max(output$Bt)/1000)* eta*
                                exp(-delta*(output$time[output$Bt==max(output$Bt)])) - 
                                c*(nstart.r["P"])  #Net profit
    lhcpars.aqua.fill[i,6] = (p*(max(output$Bt)/1000)*eta*
                                 exp(-delta*(output$time[output$Bt==max(output$Bt)])) - 
                                 c*(nstart.r["P"])) / (c*(nstart.r["P"])) # ROI 
    if(i %% 100==0) print(i)
  }
  
  hist(lhcpars.aqua.fill$h.t, breaks = 30)
  
#Next run epi mod to equilibrium over all parameter sets ##########  
  lhcpars.epi = lhcpars[,which(colnames(lhcpars) %in% names(par.snails))]
  lhcpars.epi.fill = data.frame(W = 0,
                                S1 = 0, S2 = 0, S3 = 0,
                                E1 = 0, E2 = 0, E3 = 0,
                                I1 = 0, I2 = 0, I3 = 0,
                                I.t = 0, N.t = 0)
  
  for(i in 1:nsims){
    parsim = lhcpars.epi[i,]
    
    #Run epi mod to eqbm by itself, store key outcomes
    sn.run = ode(nstart.sn,t.sn,snail_epi_allvh,parsim)
    sn.eqbm = sn.run[dim(sn.run)[1],]
    
    lhcpars.epi.fill[i,1] = sn.eqbm['Wt']
    lhcpars.epi.fill[i,2:10] = sn.eqbm[2:10]
    
    if(i %% 100==0) print(i)
    
  }
  
  lhcpars.epi.fill$I.t = rowSums(lhcpars.epi.fill[,8:10])
  lhcpars.epi.fill$N.t = rowSums(lhcpars.epi.fill[,2:10])
  
  
#Next use epi mod equilibrium values, get stocking events based on aquaculture mod run above, then run intervention over 10 yrs ##########  
  lhcpars.all.fill = data.frame(W.all = rep(0,1000),
                                I.t.all = rep(0,1000),
                                N.t.all = rep(0,1000))
  
  lhcpars.all.sims = array(data = NA, dim = c(3651, 14, nsims))
  
clus1 = makeCluster(detectCores()-1)
  clusterExport(clus1, c('lhcpars.all.fill', 'lhcpars.all.sims', 'lhcpars.aqua', 'lhcpars.epi', 'lhcpars.epi.fill', 'lhcpars',
                         'snail_epi_allvh', 'prawn_biomass'))
  
  
  for(i in 1:nsims){
    sn.eqbm = lhcpars.epi.fill[i,]
    #create new starting conditions based on the parameter set 
    nstart.lhc = unlist(c(sn.eqbm['S1'], sn.eqbm['S2'], sn.eqbm['S3'], 
                   sn.eqbm['E1'], sn.eqbm['E2'], sn.eqbm['E3'], 
                   sn.eqbm['I1'], sn.eqbm['I2'], sn.eqbm['I3'], 
                   sn.eqbm['W'], sn.eqbm['W'],
                   P = opt.ros$P_nought, L = opt.ros$L_nought)) #Optimal stocking parameters from rosenbergii ROI optimization
    
    #draw stocking parameters from prior aquaculture runs
    harvest.t.lhc = lhcpars.aqua.fill[i,1]
    n.harvest.lhc = floor(max(t.all)/harvest.t.lhc)
    stocks.lhc = data.frame(var = rep(c('P', 'L'), n.harvest.lhc),
                            time = rep(harvest.t.lhc*c(1:(n.harvest.lhc)), each = 2),
                            value = rep(c(nstart.lhc['P'], nstart.lhc['L']), n.harvest.lhc),
                            method = rep('rep', n.harvest.lhc*2))
    
    #get full parameter set from latin hypercube
    par.all.lhc = lhcpars[i,-which(colnames(lhcpars) %in% c('c', 'p', 'eta'))]
    
    #simulate
    sim.lhc = as.data.frame(ode(nstart.lhc,t.all,snail_prawn_model,par.all.lhc,
                                events = list(data = stocks.lhc)))
    
    lhcpars.all.fill[i,1] = median(sim.lhc$W[c(((years-1)*365):(years*365))+1])
    
    lhcpars.all.fill[i,2] = median(sim.lhc$I1[c(((years-1)*365):(years*365))+1],
                                   sim.lhc$I2[c(((years-1)*365):(years*365))+1] + 
                                   sim.lhc$I3[c(((years-1)*365):(years*365))+1])
    
    lhcpars.all.fill[i,3] = median(sim.lhc$I1[c(((years-1)*365):(years*365))+1],
                                   sim.lhc$I2[c(((years-1)*365):(years*365))+1] + 
                                   sim.lhc$I3[c(((years-1)*365):(years*365))+1] +
                                   sim.lhc$E3[c(((years-1)*365):(years*365))+1] +
                                   sim.lhc$E2[c(((years-1)*365):(years*365))+1] +
                                   sim.lhc$E1[c(((years-1)*365):(years*365))+1] +
                                   sim.lhc$S3[c(((years-1)*365):(years*365))+1] +
                                   sim.lhc$S2[c(((years-1)*365):(years*365))+1] +
                                   sim.lhc$S1[c(((years-1)*365):(years*365))+1])
    
    lhcpars.all.sims[ , , i] = round(as.matrix(sim.lhc), digits = 2)
    
     if(i %% 50==0) print(c(i, lhcpars.all.fill[i,]))
    
  }
  
stopCluster(clus1)

lhcfin = cbind(lhcpars, lhcpars.aqua.fill,lhcpars.epi.fill, lhcpars.all.fill)

save.image("~/RemaisWork/Schisto/Stanford/Prawn_fisheries_Public_health/Sensitivity_Analysis/sens_sims_no_diag.RData")