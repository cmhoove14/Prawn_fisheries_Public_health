#Assess sensitivity of each outcome ############
#  source('Combined_Model/epi_no_diag_prawn_mod.R')
  source('Prawn_aquaculture/prawn_aquaculture_mod.R')
#  source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
  source('Epi_Model/snail_epi_mod_no_diag_immigration.R')
  source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')

#  load('Prawn_aquaculture/aquaculture_sims.RData')
  
#get all parameters in the same vector ###########
par.all = c(par.aqua, 
            c = cost, 
            p = price, 
            eta = as.numeric(predict(eta.lm, newdata = data.frame(dens = 2.6/area))), 
            delta = -log(1-0.03)/365,
            par.snails.imm,
            eps = sqrt(area),   # Attack rate penalty 
            a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  
#generate parameter ranges to sample over ################# 
  nsims = 1000
  
  paranges = matrix(ncol = length(par.all), nrow = nsims)
  paranges = as.data.frame(paranges)
  colnames(paranges) = names(par.all)
  
  # Fill aquaculture-relevant parameter ranges
    paranges$a.p = seq(-2.7180, -2.5076, length.out = nsims)
    paranges$b.p = seq(3.4617, 3.6386, length.out = nsims)
    paranges$gam = seq(6.5e-6, 3.53e-6, length.out = nsims)
    paranges$muP = seq(0.0054, 0.0122, length.out = nsims)
    paranges$d = seq(-0.461, -0.289, length.out = nsims)
    paranges$om = seq(7.5e-9, 3.5e-9, length.out = nsims)
    paranges$k = seq(1.24/365, 3.19/365, length.out = nsims)
    paranges$linf = seq(184, 234, length.out = nsims)
    paranges$k.ros = seq(0.235/30, 0.371/30, length.out = nsims)
    paranges$c = seq(0.045, 0.12, length.out = nsims)
    paranges$p = seq(11, 22, length.out = nsims)
    paranges$eta = seq(0.24, 0.58, length.out = nsims)
    paranges$delta = seq(-log(1-0.01)/365, -log(1-0.05)/365, length.out = nsims)
    
    paranges$H = rep(1.5*area, nsims)     #don't vary these. They're only included as parameters for modeling convenience, but we know what they are
    paranges$A = rep(area, nsims)
    
  for(i in 16: length(par.all)){
    paranges[,i] = seq(par.all[i]*0.75, par.all[i]*1.25, length.out = nsims)
  }
  
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
  