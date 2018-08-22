#Assess sensitivity of each outcome ############
  source('Combined_Model/epi_no_diag_prawn_mod.R')
  source('Prawn_aquaculture/prawn_aquaculture_mod.R')
  source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
  source('Epi_Model/snail_epi_mod_no_diag.R')
  load('Prawn_aquaculture/aquaculture_sims.RData')
  
#get all parameters in the same vector ###########
par.all = c(par.aqua, par.snails,
            a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014

  par.all['c'] = cost
  par.all['p'] = p
  par.all['eta'] = eta.ros

  
#generate parameter ranges to sample over ################# 
  nsims = 1000
  
  paranges = matrix(ncol = length(par.all), nrow = nsims)
  paranges = as.data.frame(paranges)
  colnames(paranges) = names(par.all)
  
  for(i in 1: length(par.all)){
    paranges[,i] = seq(par.all[i]*0.75, par.all[i]*1.25, length.out = nsims)
  }
  
  #don't vary these. They're only included as parameters for modeling convenience, but we know what they are
  paranges$H = rep(1500, nsims)
  paranges$A = rep(1000, nsims)
  
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
  