#Assess sensitivity of each outcome ############
source('Combined_Model/epi_no_diag_prawn_mod.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
source('Epi_Model/snail_epi_mod_no_diag.R')
load('Prawn_aquaculture/aquaculture_sims.RData')

library(sensitivity)
library(parallel)

nsims = 1000
years = 10
t.all = c(0:(years*365))
nstart.r = c(P = opt.ros$P_nought, L = opt.ros$L_nought) 
eta = predict(eta.lm, newdata = data.frame(dens = nstart.r['P']/area))

#get all parameters in the same vector ###########
par.all = c(par.aqua, par.snails,
            a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014

  par.all['c'] = 0.1
  par.all['p'] = 12
  par.all['eta'] = eta
 
#generate parameter ranges to sample over ################# 
  paranges = matrix(ncol = length(par.all), nrow = nsims)
  paranges = as.data.frame(paranges)
  colnames(paranges) = names(par.all)
  
  for(i in 1: length(par.all)){
    paranges[,i] = seq(par.all[i]*0.75, par.all[i]*1.25, length.out = nsims)
  }
  
  #don't vary these. They're only included as parameters for modeling convenience, but we know what they are
  paranges$H = rep(1500, nsims)
  paranges$A = rep(1000, nsims)
  
#Check monotinicity of outcomes wrt to parameter ranges ################
  #only interested in Epi parameters, reduce dimensionality to check smaller range as we're more interested in trends
  mono_pars_epi <- paranges[round(seq(1,1000, length.out = 21)), which(colnames(paranges) %in% names(par.snails))]
  
  mono_fill_epi <- array(data = NA, dim = c(nrow(mono_pars_epi), dim(allvh.eqbm)[2], ncol(mono_pars_epi)-2))
  
  #get_eqbm <- function(par){
  for(j in 3:ncol(mono_pars_epi)){ #skip Area and human density, start at column 3 rather than column 2
      for(i in 1:nrow(mono_pars_epi)){
        test_par = colnames(mono_pars_epi)[j]
        pars_use = c(mono_pars_epi[i , which(colnames(mono_pars_epi) == test_par)], par.snails[which(names(par.snails) != test_par)])
        names(pars_use)[1] <- test_par
        sn.run = ode(nstart.sn, t.sn, snail_epi_allvh, pars_use)
        mono_fill_epi[i, 1:10,j-2] = sn.run[dim(sn.run)[1],2:11]
      } 
    
    print(j)
    
  }  
  
  mono_fill_epi_copy<-mono_fill_epi
  
#Add infected snail density and total snail density  
  for(i in 1:ncol(mono_pars_epi)-2){
    mono_fill_epi[ , 11 ,i] <- rowSums(mono_fill_epi[ , 1:9 ,i])
    mono_fill_epi[ , 12 ,i] <- rowSums(mono_fill_epi[ , 7:9 ,i])
  }

#Plots to check monotinicity    
  par(mfrow = c(3,1), mar = c(4,3,1,1))
  
  for(i in 1:(ncol(mono_pars_epi)-2)){
    plot(mono_pars_epi[, (i+2)], mono_fill_epi[ , 10, i], xlab = colnames(mono_pars_epi[(i+2)]), ylab = "W",
         type = "l", col = "purple", lwd = 2)
    plot(mono_pars_epi[, (i+2)], mono_fill_epi[ , 11, i], xlab = colnames(mono_pars_epi[(i+2)]), ylab = "N",
         type = "l", col = "Black", lwd = 2)
    plot(mono_pars_epi[, (i+2)], mono_fill_epi[ , 12, i], xlab = colnames(mono_pars_epi[(i+2)]), ylab = "I",
         type = "l", col = "Red", lwd = 2)
    print(i)
  }
  
  par(mfrow = c(1,1))
  
#  Relationship between muN2 and total snail density is a bit concave, but otherwise everything checks out as monotonic
  
#Augment latin hypercube to sample from ###########
  #create a matrix of indices for the LHC where nrow=number of sims and ncol=number of variables
  LHC_indices<-matrix(0, nrow=nsims, ncol=length(par.all))  
  #Add structure of latin hypercube
  set.seed(43018)
  for(j in 1:dim(LHC_indices)[2]){
    LHC_indices[,j]<-sample(1:nsims, size=nsims, replace=FALSE)
  }  
  
  lhcpars = paranges
  
#reorder parameters in order of LHC indices  
  for(k in 1:length(par.all)){
    lhcpars[,k] = paranges[,k][LHC_indices[,k]]
  }
  
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