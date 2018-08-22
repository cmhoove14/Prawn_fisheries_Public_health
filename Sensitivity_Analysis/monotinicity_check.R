source('Combined_Model/epi_no_diag_prawn_mod.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
source('Epi_Model/snail_epi_mod_no_diag.R')
source('Sensitivity_Analysis/sensitivity_prep.R')
load('Prawn_aquaculture/aquaculture_sims.RData')


#Model prep
years = 10
t.all = c(0:(years*365))
nstart.r = c(P = opt.ros$P_nought, L = opt.ros$L_nought) 
eta = predict(eta.lm, newdata = data.frame(dens = nstart.r['P']/area))

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
