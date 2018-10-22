#Assess sensitivity of each outcome ############
#source("Data/aquaculture_parameters.R")
#source("Data/snail_epi_parameters.R")
#source("Data/combined_parameters.R")

#get all parameters in the same vector ###########
par.all = c(par.aqua, 
            c = cost, 
            fc = fixed_cost,
            p = price, 
            eta = as.numeric(predict(eta.lm, newdata = data.frame(dens = 2.6/area))), 
            delta = -log(1-0.07)/365,
            par.snails.imm, 
            par.epi_prawn,
            weight_hi = 0.05,
            weight_lo = 0.014,
            eff = 0.85,
            epmL = 3.6)

#generate parameter ranges to sample over ################# 
  nsims = 1000
  
  paranges = as.data.frame(matrix(ncol = length(par.all), nrow = nsims))
    colnames(paranges) = names(par.all)
  
  # Fill aquaculture-relevant parameter ranges
    paranges$a.p = seq(-2.7180, -2.5076, length.out = nsims)
    paranges$b.p = seq(3.4617, 3.6386, length.out = nsims)
    paranges$gam = seq(6.5e-6, 3.5e-6, length.out = nsims)
    paranges$muP = seq(2/365, 3/365, length.out = nsims)
    paranges$d = seq(-0.461, -0.289, length.out = nsims)
    paranges$om = seq(7.5e-9, 3.5e-9, length.out = nsims)
    paranges$k = seq(1.24/365, 3.19/365, length.out = nsims)
    paranges$linf = seq(184, 234, length.out = nsims)
    paranges$k.ros = seq(0.235/30, 0.371/30, length.out = nsims)
    paranges$c = seq(0.045, 0.12, length.out = nsims)
    paranges$fc = 0 # seq(0, 1000, length.out = nsims) influence of fixed cost explored in separate analysis
    paranges$p = seq(11, 22, length.out = nsims)
    paranges$eta = seq(0.24, 0.58, length.out = nsims)
    paranges$delta = seq(-log(1-0.03)/365, -log(1-0.20)/365, length.out = nsims)
    
    paranges$H = rep(1.5*area, nsims)     #don't vary these. They're only included as parameters for modeling convenience, but we know what they are
    paranges$A = rep(area, nsims)

  # Fill epi-relevant parameter ranges
    paranges$f = seq(0.06, 0.36, length.out = nsims)
    paranges$Kn = seq(25*area, 75*area, length.out = nsims)
    paranges$z = seq(0.25, 1, length.out = nsims)
    paranges$xi = seq(0, 1/365, length.out = nsims)
    paranges$muN1 = seq(1/25, 1/100, length.out = nsims)
    paranges$muN2 = seq(1/50, 1/125, length.out = nsims)
    paranges$muN3 = seq(1/75, 1/150, length.out = nsims)
    paranges$muI = seq(1/7, 1/20, length.out = nsims)
    paranges$nfr = seq(1.1, 4, length.out = nsims)
    paranges$g1 = seq(1/20, 1/60, length.out = nsims)
    paranges$g2 = seq(1/40, 1/100, length.out = nsims)
    paranges$beta = seq(2e-7, 8e-7, length.out = nsims)
    paranges$m = seq(0.5, 1.5, length.out = nsims)
    paranges$sigma = seq(1/70, 1/30, length.out = nsims)
    paranges$lambda = seq(5e-6, 1e-5, length.out = nsims)
    paranges$theta1 = seq(1, 5, length.out = nsims)
    paranges$theta2 = seq(2, 10, length.out = nsims)
    paranges$phi = seq(0.01, 0.3, length.out = nsims)
    paranges$muW = seq(1/(2*365), 1/(4*365), length.out = nsims)
    paranges$muH = seq(1/(50*365), 1/(70*365), length.out = nsims)
    paranges$psi1 = 0
    paranges$psi2 = 0
    paranges$psi3 = 0
    paranges$siteS1 = 0
    paranges$siteS2 = 0
    paranges$siteS3 = 0
    paranges$siteE1 = 0
    paranges$siteE2 = 0
    paranges$siteE3 = 0
    paranges$siteI1 = 0
    paranges$siteI2 = 0
    paranges$siteI3 = 0
    
  #Fill combined model parameter ranges
    paranges$eps = seq(1, 100, length.out = nsims)
    paranges$a.s = seq(0.1, 0.3, length.out = nsims)
    paranges$b.s = seq(2, 3, length.out = nsims)
    paranges$ar.slp = seq(0.5, 1.5, length.out = nsims)
    paranges$th = seq(0.2, 0.5, length.out = nsims)
    
  #fill treatement/outcomes related parameters
    paranges$weight_hi = seq(0.03, 0.15, length.out = nsims)
    paranges$weight_lo = seq(0.003,0.03, length.out = nsims)
    paranges$cvrg = seq(0.5, 0.95, length.out = nsims)
    paranges$eff = seq(0.75, 0.95, length.out = nsims)
    paranges$epmL = seq(2, 5, length.out = nsims)
    
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
  