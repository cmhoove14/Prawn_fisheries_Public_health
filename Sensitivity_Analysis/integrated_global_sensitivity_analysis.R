#Assess sensitivity of each outcome ############
source('Combined_Model/epi_prawn_mod.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
source('Epi_Model/snail_epi_mod.R')

library(sensitivity)
library(parallel)

nsims = 1000
years = 10
t.all = c(0:(years*365))
nstart.r = c(P = 19000, L = 38) 
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
    paranges[,i] = seq(par.all[i]*0.7, par.all[i]*1.3, length.out = nsims)
  }
  
  #don't vary these. They're only included as parameters for modeling convenience, but we know what they are
  paranges$H = rep(1000, nsims)
  paranges$A = rep(10000, nsims)
  
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
    parun<-lhcpars.aqua[i,1:8]
    c = lhcpars.aqua[i,9]
    p = lhcpars.aqua[i,10]
    eta = lhcpars.aqua[i,11]
      
    output = as.data.frame(ode(nstart.r, t.p, prawn_biomass, parun))
    output$B = ((par.aqua['a.p']*(output$L/10)^par.aqua['b.p'])/10)       
    output$Bt = output$B*output$P      
    
    lhcpars.aqua.fill[i,1] = output$time[output$B*output$P==max(output$B*output$P)]   #harvest time
    lhcpars.aqua.fill[i,2] = max(output$B*output$P)                       #harvest total biomass
    lhcpars.aqua.fill[i,3] = max(output$B*output$P) * eta                 #harvest marketable biomass
    lhcpars.aqua.fill[i,4] = p*(max(output$B*output$P)/1000)               # Raw Profit (revenue)
    lhcpars.aqua.fill[i,5] = p*(max(output$B*output$P)/1000)*
                                exp(-delta*(output$time[output$B*output$P==max(output$B*output$P)])) - 
                                c*(nstart.r["P"])  #Net profit
    lhcpars.aqua.fill[i,6] = (p*(max(output$B*output$P)/1000)*
                                 exp(-delta*(output$time[output$B*output$P==max(output$B*output$P)])) - 
                                 c*(nstart.r["P"])) / (c*(nstart.r["P"])) # ROI 
    if(i %% 100==0) print(i)
  }
  
  hist(lhcpars.aqua.fill$h.t, breaks = 30)
  
#Next run epi mod to equilibrium over all parameter sets ##########  
  lhcpars.epi = lhcpars[,which(colnames(lhcpars) %in% names(par.snails))]
  lhcpars.epi.fill = data.frame(W = 0,
                                S1 = 0, S2 = 0, S3 = 0,
                                E1 = 0, E2 = 0, E3 = 0,
                                I2 = 0, I3 = 0,
                                I.t = 0, N.t = 0)
  
  for(i in 1:nsims){
    parsim = lhcpars.epi[i,]
    
    #Run epi mod to eqbm by itself, store key outcomes
    sn.run = ode(nstart.sn,t.sn,snail_epi,parsim)
    sn.eqbm = sn.run[dim(sn.run)[1],]
    
    lhcpars.epi.fill[i,1] = sn.eqbm['W']
    lhcpars.epi.fill[i,2:9] = sn.eqbm[2:9]
    
    if(i %% 100==0) print(i)
    
  }
  
  lhcpars.epi.fill$I.t = rowSums(lhcpars.epi.fill[,8:9])
  lhcpars.epi.fill$N.t = rowSums(lhcpars.epi.fill[,2:9])
  
  
#Next use epi mod equilibrium values, get stocking events based on aquaculture mod run above, then run intervention over 10 yrs ##########  
  lhcpars.all.fill = data.frame(W.all = rep(0,1000),
                                I.t.all = rep(0,1000),
                                N.t.all = rep(0,1000))
  
  lhcpars.all.sims = array(data = NA, dim = c(3651, 12, nsims))
  
clus1 = makeCluster(detectCores()-1)
  clusterExport(clus1, c('lhcpars.all.fill', 'lhcpars.all.sims', 'lhcpars.aqua', 'lhcpars.epi', 'lhcpars.epi.fill', 'lhcpars',
                         'snail_epi', 'prawn_biomass'))
  
  
  for(i in 1:nsims){
    sn.eqbm = lhcpars.epi.fill[i,]
    #create new starting conditions based on the parameter set 
    nstart.lhc = c(sn.eqbm['S1'], sn.eqbm['S2'], sn.eqbm['S3'], 
                   sn.eqbm['E1'], sn.eqbm['E2'], sn.eqbm['E3'], 
                   sn.eqbm['I2'], sn.eqbm['I3'], sn.eqbm['W'], 
                   P = 19000, L = 38) #Optimal stocking parameters from rosenbergii ROI optimization
    
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
    
    lhcpars.all.fill[i,2] = median(sim.lhc$I2[c(((years-1)*365):(years*365))+1] + 
                                   sim.lhc$I3[c(((years-1)*365):(years*365))+1])
    
    lhcpars.all.fill[i,3] = median(sim.lhc$I2[c(((years-1)*365):(years*365))+1] + 
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
 #lhcfin$N.t.all.10yr = lhcpars.all.sims[3651, 10, ]
 #lhcfin$W.all.10yr = lhcpars.all.sims[3651, 10, ]
  #save(lhcfin, file='Sensitivity_Analysis/lhc_prcc_df_plusminus30%_inputs.Rdata')
  #load("~/RemaisWork/Schisto/Stanford/Prawn_fisheries_Public_health/Sensitivity_Analysis/lhc_sims_n1000_seed043093_30per_par_var.Rdata")
  #lhcpars.all.fill[,1] = lhcpars.all.sims[3651, 10, ]
  #lhcpars.all.fill[,2] = rowSums(lhcpars.all.sims[3651, 8:9, ])
  #lhcpars.all.fill[,3] = rowSums(lhcpars.all.sims[3651, 2:9, ])
# Assess PRCC of prawn model outputs ###############
# profit sensitivity  
  profit.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]),
                   y = as.data.frame(lhcpars.aqua.fill$profit),
                   rank = TRUE)
  
  profit.pcc.df = profit.pcc$PRCC
  profit.pcc.df$var = c(names(par.aqua), 'c', 'p', 'eta')
  
  profit.lhsprcc = ggplot(profit.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
    geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
    geom_hline(yintercept = 0, lty = 2)+
    labs(x = 'Parameter', y = 'Profit PRCC')
  
  profit.lhsprcc
  
# roi sensitivity  
  roi.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]),
                y = as.data.frame(lhcpars.aqua.fill$roi),
                rank = TRUE)
  
  roi.pcc.df = roi.pcc$PRCC
  roi.pcc.df$var = c(names(par.aqua), 'c', 'p', 'eta')
  
  roi.lhsprcc = ggplot(roi.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
    geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
    geom_hline(yintercept = 0, lty = 2)+
    labs(title = 'roi sensitivity', x = 'Parameter', y = 'PRCC')
  
  roi.lhsprcc
  
# harvest time sensitivity  
  ht.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% c(names(par.aqua), 'c', 'p', 'eta'))]),
               y = as.data.frame(lhcpars.aqua.fill$h.t),
               rank = TRUE)
  
  ht.pcc.df = ht.pcc$PRCC
  ht.pcc.df$var = c(names(par.aqua), 'c', 'p', 'eta')
  
  ht.lhsprcc = ggplot(ht.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
    geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(title = 'ht sensitivity', x = 'Parameter', y = 'PRCC')
  
  ht.lhsprcc
  
# Assess PRCC of epi model outputs ###########
#Equilibrium mean worm burden  
  W.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% 
                                                 names(par.snails)[-which(names(par.snails) %in% 
                                                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
              y = as.data.frame(lhcpars.epi.fill$W),
              rank = TRUE)
  
  W.pcc.df = W.pcc$PRCC
  W.pcc.df$var = names(par.snails)[-which(names(par.snails) %in% 
                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))]
  
  W.lhsprcc = ggplot(W.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
    geom_bar(fill = 'purple', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(x = 'Parameter', y = 'Equlibrium W PRCC')
  
  W.lhsprcc
  
#Equilibrium total snail population 
  Nt.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% 
                                                 names(par.snails)[-which(names(par.snails) %in% 
                                                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
              y = as.data.frame(lhcpars.epi.fill$N.t),
              rank = TRUE)
  
  Nt.pcc.df = Nt.pcc$PRCC
  Nt.pcc.df$var = names(par.snails)[-which(names(par.snails) %in% 
                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))]
  
  Nt.lhsprcc = ggplot(Nt.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
    geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(title = 'Nt sensitivity', x = 'Parameter', y = 'PRCC')
  
  Nt.lhsprcc
  
  pairs(cbind(lhcpars.epi.fill$N.t, lhcpars.epi[,3:8]))
  pairs(cbind(lhcpars.epi.fill$N.t, lhcpars.epi[,13:22]))
  
#Equilibrium infected snail population 
  It.pcc = pcc(X = as.data.frame(lhcpars[,which(colnames(lhcpars) %in% 
                                                  names(par.snails)[-which(names(par.snails) %in% 
                                                                             c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
               y = as.data.frame(lhcpars.epi.fill$I.t),
               rank = TRUE)
  
  It.pcc.df = It.pcc$PRCC
  It.pcc.df$var = names(par.snails)[-which(names(par.snails) %in% 
                                             c('A', 'H', 'psi1', 'psi2', 'psi3'))]
  
  It.lhsprcc = ggplot(It.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5,0,0.5))+
    geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(x = 'Parameter', y = 'Infected snail density PRCC')
  
  It.lhsprcc
  
# Assess PRCC of combined model (with intervention regime) outputs ###########
# Ending mean worm burden  
  Wall.pcc = pcc(X = as.data.frame(lhcfin[,which(colnames(lhcfin) %in% 
                                                 names(par.all)[-which(names(par.all) %in% 
                                                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
              y = as.data.frame(lhcfin$W.all),
              rank = TRUE)
  
  Wall.pcc.df = Wall.pcc$PRCC
  Wall.pcc.df$var = names(par.all)[-which(names(par.all) %in% 
                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))]
  
  Wall.lhsprcc = ggplot(Wall.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-0.75,0.75), breaks = c(-0.75,-0.25,-0.5,0,0.25,0.5,0.75))+
    geom_bar(fill = 'purple', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(x = 'Parameter', y = 'Final W PRCC')
  
  Wall.lhsprcc
  
#Ending infected snail density  
  I_tall.pcc = pcc(X = as.data.frame(lhcfin[,which(colnames(lhcfin) %in% 
                                                   names(par.all)[-which(names(par.all) %in% 
                                                                           c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
                 y = as.data.frame(lhcfin$I.t.all),
                 rank = TRUE)
  
  I_tall.pcc.df = I_tall.pcc$PRCC
  I_tall.pcc.df$var = names(par.all)[-which(names(par.all) %in% 
                                            c('A', 'H', 'psi1', 'psi2', 'psi3'))]
  
  I_tall.lhsprcc = ggplot(I_tall.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1))+
    geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(title = 'Infected snail density sensitivity', x = 'Parameter', y = 'PRCC')
  
  I_tall.lhsprcc
  
#Ending total snail density  
  N_tall.pcc = pcc(X = as.data.frame(lhcfin[,which(colnames(lhcfin) %in% 
                                                     names(par.all)[-which(names(par.all) %in% 
                                                                             c('A', 'H', 'psi1', 'psi2', 'psi3'))])]),
                   y = as.data.frame(colSums(lhcpars.all.sims[3651, 2:9, ])),
                   rank = TRUE)
  
  N_tall.pcc.df = N_tall.pcc$PRCC
  N_tall.pcc.df$var = names(par.all)[-which(names(par.all) %in% 
                                              c('A', 'H', 'psi1', 'psi2', 'psi3'))]
  
  N_tall.lhsprcc = ggplot(N_tall.pcc.df, aes(x = var, y = original)) +
    theme_bw()+
    scale_y_continuous(limits = c(-0.75,0.75), breaks = c(-0.75,-0.25,-0.5,0,0.25,0.5,0.75))+
    geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
    #geom_errorbar(x = var, )
    geom_hline(yintercept = 0, lty = 2)+
    labs(x = 'Parameter', y = 'Final snail density PRCC')
  
  N_tall.lhsprcc
  
#Something weird going on here, PRCC is near 0 for all variables
  plot(lhcfin$N.t, lhcfin$N.t.all, pch = 16, cex = 0.7)
  plot(lhcfin$N.t.all, colSums(lhcpars.all.sims[3651, 2:9, ]), pch = 16, cex = 0.7,
       xlab = 'fin snail pop', ylab = 'median last year snail pop')
  
  lhcfin$N.t.all.10yr = colSums(lhcpars.all.sims[3651, 2:9, ])
    pairs(lhcfin[,c(1:8, 49)], pch = 18, cex = 0.7)
    pairs(lhcfin[,c(1:8, 50)], pch = 18, cex = 0.7)
    
    pairs(lhcfin[,c(11:17, 50)], pch = 18, cex = 0.7)
    pairs(lhcfin[,c(21:34, 50)], pch = 18, cex = 0.7)
    
# Combine plots of key model outcomes ##################
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    library(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                        layout.pos.col = matchidx$col))
      }
    }
  }
  
  fig1.layout = matrix(c(1,2,
                         3,3), ncol = 2, byrow = T)
  
  windows(width = 300, height = 200)
  multiplot(profit.lhsprcc, It.lhsprcc, 
            Wall.lhsprcc, 
            layout = fig1.layout)
  
##################
##################