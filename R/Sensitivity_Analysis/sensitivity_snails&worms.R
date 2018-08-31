#Assess sensitivity of each outcome ############
source('Combined_Model/epi_prawn_mod.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
source('Epi_Model/snail_epi_mod.R')

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
    paranges[,i] = seq(par.all[i]*0.75, par.all[i]*1.25, length.out = nsims)
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
  
#Next run epi mod to equilibrium, get stocking events based on aquaculture mod run above, then run intervention over 10 yrs ##########  
  lhcpars.epi = lhcpars[,which(colnames(lhcpars) %in% names(par.snails))]
  lhcpars.epi.fill = data.frame(W = 0,
                                I.t = 0,
                                N.t = 0)
  
  lhcpars.all.fill = data.frame(W = 0,
                                I.t = 0,
                                N.t = 0)
  
  for(i in 1:nsims){
    #Run epi mod to eqbm by itself, store key outcomes
    sn.run = as.data.frame(ode(nstart.sn,t.sn,snail_epi,par.snails))
    sn.eqbm = sn.run[dim(sn.run)[1],]
    
    lhcpars.epi.fill[i,1] = sn.eqbm$W
    lhcpars.epi.fill[i,2] = sn.eqbm$I2 + sn.eqbm$I3
    lhcpars.epi.fill[i,3] = sum(sn.eqbm[,2:9])
    
    #create new starting conditions based on the parameter set 
    nstart.lhc = c(S1 = sn.eqbm$S1, S2 = sn.eqbm$S2, S3 = sn.eqbm$S3, 
                   E1 = sn.eqbm$E1, E2 = sn.eqbm$E2, E3 = sn.eqbm$E3, 
                   I2 = sn.eqbm$I2, I3 = sn.eqbm$I3, W = sn.eqbm$W, 
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
    print(i)
    # if(i %% 100==0) print(i)
    
  }
  
lhcfin = cbind(lhcpars, lhcpars.aqua.fill,lhcpars.epi.fill, lhcpars.all.fill)
  colnames(lhcfin)[47:49] = c('W.all', 'I.t.all', 'N.t.all')
  save(lhcfin, file='Sensitivity_Analysis/lhc_prcc_df.Rdata')
#check model behavior ###############
  harvest.time = harvest.time
  ncycles = 10
  nstart.lt = nstart
  nstart.lt['P'] = 0
  nstart.lt['L'] = 0
  nstart.lt['W'] = 72
  
  time.lt = seq(0, harvest.time*ncycles, 1)

# Define stocking events at the beginning of each cycle
stocking = data.frame(var = rep(c('P', 'L'), ncycles),
                      time = rep(seq(pzq.delay, pzq.delay + harvest.time*(ncycles-1), harvest.time), each = 2),
                      value = rep(c(nstart['P'], nstart['L']), ncycles),
                      method = rep('rep', ncycles*2))

output.lt = as.data.frame(ode(nstart.lt, time.lt, snail_prawn_model, parameters,
                              events = list(data = stocking)))
#Store end values to use as start below  
eqbm = output.lt[dim(output.lt)[2],c(2:ncol(output.lt))]

  output.lt$S.t = output.lt$S1 + output.lt$S2 + output.lt$S3                      # Total susceptible snails
  output.lt$E.t = output.lt$E1 + output.lt$E2 + output.lt$E3                      # Total exposed snails 
  output.lt$I.t = output.lt$I2 + output.lt$I3                                     # Total infected snails
  output.lt$N1.t = output.lt$S1 + output.lt$E1                                    # Total snails of size class 1
  output.lt$N2.t = output.lt$S2 + output.lt$E2 + output.lt$I2                     # Total snails of size class 2
  output.lt$N3.t = output.lt$S3 + output.lt$E3 + output.lt$I3                     # Total snails of size class 3
  output.lt$N.t = output.lt$S.t + output.lt$E.t + output.lt$I.t                   # Total snails
  output.lt$prev = pnbinom(2, size = 0.2, mu = output.lt$W, lower.tail = FALSE)   # Estimated prevalence, using a negative binomial dist. with k = 0.2 (fitted from EPLS data)
  output.lt$B = ((parameters['a.p']*(output.lt$L/10)^parameters['b.p'])/10)       # Mean prawn biomass, transformed from length using allometric equation
  output.lt$Bt = output.lt$B*output.lt$P                                          # Total prawn biomass

# Plot snail dynamics by size class
plot(x = output.lt$time, y = output.lt$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output.lt$N.t)),
     main = 'Snail Size Classes')
  lines(output.lt$time, output.lt$N1.t, col = 'green', lwd = 2)
  lines(output.lt$time, output.lt$N2.t, col = 'blue', lwd = 2)
  lines(output.lt$time, output.lt$N3.t, col = 'red', lwd = 2)
  lines(output.lt$time, output.lt$P, col = 'gold2', lty=2, lwd=2)
  legend('topright', legend = c('total', '1', '2', '3'), lwd = 2, col = c('black', 'green', 'blue', 'red'), cex = 0.7)

# Plot snail dynamics by infection class
plot(x = output.lt$time, y = output.lt$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output.lt$N.t)), 
     main = 'Snail Infection Classes')
  lines(output.lt$time, output.lt$I.t, col = 'red', lwd = 2)
  lines(output.lt$time, output.lt$E.t, col = 'orange', lwd = 2)
  lines(output.lt$time, output.lt$S.t, col = 'green', lwd = 2)
  legend('topright', legend = c('total', 'S', 'E', 'I'), lwd = 2, col = c('black', 'green', 'orange', 'red'), cex = 0.7)

# Plot mean human worm burden and prevalence of infection
plot(x = output.lt$time, y = output.lt$W, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Worm burden', ylim = c(0,max (output.lt$W)),
     main = 'Worm Burden')
plot(x = output.lt$time, y = output.lt$prev, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Prevalence', ylim = c(0,max (output.lt$prev)),
     main = 'Estimated Prevalence')

#develop parameter sets in LHC to sample ############
  sims2 = 150
#snail parameters
  f.range<-seq(0.05, 0.75, length.out = sims2)
  Kn.range<-seq(10, 60, length.out = sims2)
  z.range<-seq(0, 1, length.out = sims2)
  a.s.range<-seq(parameters['a.s']/2, parameters['a.s']*2, length.out = sims2)
  b.s.range<-seq(parameters['b.s']/2, parameters['b.s']*2, length.out = sims2)
  g1.range<-seq(parameters['g1']/2, parameters['g1']*2, length.out = sims2)
  g2.range<-seq(parameters['g2']/2, parameters['g2']*2, length.out = sims2)
  muN1.range<-seq(parameters['muN1']/1.2, parameters['muN1']*1.2, length.out = sims2)
  muN2.range<-seq(parameters['muN2']/1.2, parameters['muN2']*1.2, length.out = sims2)
  muN3.range<-seq(parameters['muN3']/1.2, parameters['muN3']*1.2, length.out = sims2)
  muI.range<-seq(parameters['muI']/1.2, parameters['muI']*1.2, length.out = sims2)
  s.range<-seq(parameters['s']/10, parameters['s']*10, length.out = sims2)
  beta.range<-seq(parameters['beta']/100, parameters['beta']*100, length.out = sims2)
  lambda.range<-seq(parameters['lambda']/100, parameters['lambda']*10, length.out = sims2)
  sigma.range<-seq(1/15, 1/115, length.out = sims2)
  theta.range<-seq(2, 10, length.out = sims2)
  m.range<-seq(parameters['m']/10, parameters['m']*10, length.out = sims2)

#create a matrix of indices for the LHC where nrow=number of sims and ncol=number of variables
  LHC_indices2<-matrix(0, nrow=sims2, ncol=17)  
#Add structure of latin hypercube
  for(j in 1:dim(LHC_indices2)[2]){
    LHC_indices2[,j]<-sample(1:sims2, size=sims2, replace=FALSE)
  }  
#Fill LHC values
params02<-cbind(f=f.range[LHC_indices2[,1]], 
                Kn=Kn.range[LHC_indices2[,2]], 
                z=z.range[LHC_indices2[,3]], 
                a.s=a.s.range[LHC_indices2[,4]], 
                b.s=b.s.range[LHC_indices2[,5]], 
                g1=g1.range[LHC_indices2[,6]], 
                g2=g2.range[LHC_indices2[,7]], 
                muN1=muN1.range[LHC_indices2[,8]],
                muN2=muN2.range[LHC_indices2[,9]],
                muN3=muN3.range[LHC_indices2[,10]],
                muI=muI.range[LHC_indices2[,11]],
                s=s.range[LHC_indices2[,12]],
                beta=beta.range[LHC_indices2[,13]],
                lambda=lambda.range[LHC_indices2[,14]],
                sigma=sigma.range[LHC_indices2[,15]],
                theta=theta.range[LHC_indices2[,16]],
                m=m.range[LHC_indices2[,17]])

#Hold all other parameters constant  
constantparams<-matrix(ncol = length(parameters), nrow = sims2)

for(i in 1:length(parameters)){
  constantparams[,i] = rep(parameters[i],sims2)
}

colnames(constantparams)<-names(parameters)
vars2<-colnames(params02)
params2<-constantparams[, -which(colnames(constantparams) %in% vars2)]

#Add in sampled values of parameters of interest from above
  params.fin2<-cbind(params2, params02)  

#Start without considering any PZQ; assess sensitivity based on prcc throughout a prawn cycle #################
  ncycles2 = 5
  
  time2 = seq(0, harvest.time*ncycles2, 1)
  
  harvest.time.red<-seq(0, harvest.time, by=11)

#Outputs of interest
  snails<-numeric()
  infs<-numeric()
  ws<-numeric()
  prev<-numeric()
#matrices to fill with prcc
  snails.fill<-matrix(nrow = ncol(LHC_indices2), ncol = length(harvest.time.red))
  infs.fill<-matrix(nrow = ncol(LHC_indices2), ncol = length(harvest.time.red))
  ws.fill<-matrix(nrow = ncol(LHC_indices2), ncol = length(harvest.time.red))
  prev.fill<-matrix(nrow = ncol(LHC_indices2), ncol = length(harvest.time.red))

  nstart2 = c(S1 = eqbm[,1],
              S2 = eqbm[,2],
              S3 = eqbm[,3],
              E1 = eqbm[,4],
              E2 = eqbm[,5],
              E3 = eqbm[,6],
              I2 = eqbm[,7],
              I3 = eqbm[,8],
              W = eqbm[,9],
              P = 0,
              L = 0)


#Run model 
for(j in 1:length(harvest.time.red)){
  for(i in 1:sims2){
    print(c(i,j))
    
    parameters2<-params.fin2[i,]
    
    output2 = as.data.frame(ode(nstart2, time2, snail_prawn_model, parameters2,
                                events = list(data = rbind(stocking, pzq))))
    
    infs[i] = output2$I2[dim(output2)[1]-harvest.time.red[j]] + output2$I3[dim(output2)[1]-harvest.time.red[j]] # Total infected snails
    
    snails[i] = (output2$S1[dim(output2)[1]-harvest.time.red[j]] + 
                   output2$S2[dim(output2)[1]-harvest.time.red[j]] +
                   output2$S3[dim(output2)[1]-harvest.time.red[j]] +
                   output2$E1[dim(output2)[1]-harvest.time.red[j]] +
                   output2$E2[dim(output2)[1]-harvest.time.red[j]] +
                   output2$E3[dim(output2)[1]-harvest.time.red[j]] +
                   output2$I2[dim(output2)[1]-harvest.time.red[j]] +
                   output2$I3[dim(output2)[1]-harvest.time.red[j]])                    # Total snails
    
    ws[i] = output2$W[dim(output2)[1]-harvest.time.red[j]]                           # Mean worm burden
    
    prev[i] = pnbinom(2, size = 0.2, 
                      mu = output2$W[dim(output2)[1]-harvest.time.red[j]], 
                      lower.tail = FALSE)     # Estimated prevalence
    
    par(mfrow=c(2,1), mar = c(4,3.75,1,0.4)+0.1)  
    plot(x = output2$time, y = (output2$S1 + output2$S2 + output2$S3 +
                                  output2$E1 + output2$E2 + output2$E3 +
                                  output2$I2 + output2$I3), 
         type = 'l', col = 'black', lwd=2, ylab = 'Snails', xlab = 'time')
    lines(x = output2$time, y = output2$I2, col = 'red', lwd=2)
    
    plot(output2$time, output2$P/100, type = 'l', lwd = 2, col = 'gold2', 
         xlab = '', ylab = 'Prawn & worm pop') 
    lines(x = output2$time, y = output2$W, col = 'purple', lwd=2)
    
    #add outcome parameters of interest to corresponding parameter sets in data frame
    params.snail<-as.data.frame(cbind(params.fin2, infs, snails, ws, prev))
    ph.test2<-cbind(params.snail[, -which(colnames(params.snail) %in% colnames(params2))]) 
    
    prev.pcc<-pcc(X = as.data.frame(ph.test2[,c(1:length(vars2))]),
                  y = as.data.frame(ph.test2[,c(dim(ph.test2)[2])]),
                  rank = TRUE)
    
    prev.fill[,j]<-prev.pcc$PRCC[[1]]
    
    ws.pcc<-pcc(X = as.data.frame(ph.test2[,c(1:length(vars2))]),
                y = as.data.frame(ph.test2[,c(dim(ph.test2)[2]-1)]),
                rank = TRUE)
    
    ws.fill[,j]<-ws.pcc$PRCC[[1]]
    
    snails.pcc<-pcc(X = as.data.frame(ph.test2[,c(1:length(vars2))]),
                    y = as.data.frame(ph.test2[,c(dim(ph.test2)[2]-2)]),
                    rank = TRUE)
    
    snails.fill[,j]<-snails.pcc$PRCC[[1]]
    
    infs.pcc<-pcc(X = as.data.frame(ph.test2[,c(1:length(vars2))]),
                  y = as.data.frame(ph.test2[,c(dim(ph.test2)[2]-3)]),
                  rank = TRUE)
    
    infs.fill[,j]<-infs.pcc$PRCC[[1]]
  }
  
}  

#Analyze results to see which outcomes are sensitive to which variables; if sensitivity changes over prawn cycle    
#Process data frames ###########
snails.fill<-as.data.frame(snails.fill)  
colnames(snails.fill)<-harvest.time.red
snails.fill$mean<-numeric(length=17)
snails.fill$sd<-numeric(length=17)

for(i in 1:nrow(snails.fill)){
  snails.fill[i, length(harvest.time.red)+1] = mean(as.numeric(snails.fill[i, c(1:length(harvest.time.red))]))
  snails.fill[i, length(harvest.time.red)+2] = sd(as.numeric(snails.fill[i, c(1:length(harvest.time.red))]))
}
snails.fill$vars<-vars2

infs.fill<-as.data.frame(infs.fill)  
colnames(infs.fill)<-harvest.time.red
infs.fill$mean<-numeric(length=17)
infs.fill$sd<-numeric(length=17)

for(i in 1:nrow(infs.fill)){
  infs.fill[i, length(harvest.time.red)+1] = mean(as.numeric(infs.fill[i, c(1:length(harvest.time.red))]))
  infs.fill[i, length(harvest.time.red)+2] = sd(as.numeric(infs.fill[i, c(1:length(harvest.time.red))]))
}
infs.fill$vars<-vars2  

ws.fill<-as.data.frame(ws.fill)  
colnames(ws.fill)<-harvest.time.red
ws.fill$mean<-numeric(length=17)
ws.fill$sd<-numeric(length=17)

for(i in 1:nrow(ws.fill)){
  ws.fill[i, length(harvest.time.red)+1] = mean(as.numeric(ws.fill[i, c(1:length(harvest.time.red))]))
  ws.fill[i, length(harvest.time.red)+2] = sd(as.numeric(ws.fill[i, c(1:length(harvest.time.red))]))
}
ws.fill$vars<-vars2    

prev.fill<-as.data.frame(prev.fill)  
colnames(prev.fill)<-harvest.time.red
prev.fill$mean<-numeric(length=17)
prev.fill$sd<-numeric(length=17)

for(i in 1:nrow(prev.fill)){
  prev.fill[i, length(harvest.time.red)+1] = mean(as.numeric(prev.fill[i, c(1:length(harvest.time.red))]))
  prev.fill[i, length(harvest.time.red)+2] = sd(as.numeric(prev.fill[i, c(1:length(harvest.time.red))]))
}
prev.fill$vars<-vars2  

#Plot mean +/- sd prcc of each outcome ###############     
ggplot(snails.fill, aes(x = vars, y = mean)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_bar(fill = 'blue', stat = 'identity', width = 0.25)+
  labs(title = 'Snail population sensitivity', x = 'Parameter', y = 'Mean PRCC')+
  geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)))

ggplot(infs.fill, aes(x = vars, y = mean)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_bar(fill = 'red', stat = 'identity', width = 0.25)+
  labs(title = 'Infected snails sensitivity', x = 'Parameter', y = 'Mean PRCC')+
  geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.1)

ggplot(ws.fill, aes(x = vars, y = mean)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_bar(fill = 'purple', stat = 'identity', width = 0.25)+
  labs(title = 'Mean worm burden sensitivity', x = 'Parameter', y = 'Mean PRCC')+
  geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.1)

ggplot(prev.fill, aes(x = vars, y = mean)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
  labs(title = 'Prevalence sensitivity', x = 'Parameter', y = 'Mean PRCC')+
  geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.1)

#Plot prcc estimates across the prawn cycle to see if/how outcome sensitivity changes ############
plot(x = as.numeric(colnames(snails.fill)), y = snails.fill[1,c(1:33)], type = 'l',
     ylim = c(-1,1), xlab = 'Harvest cycle', ylab = 'prcc')
for(i in 2:nrow(snails.fill)){
  lines(x = as.numeric(colnames(snails.fill)[1:33]), y = snails.fill[i,c(1:33)], col = i)
}

plot(x = as.numeric(colnames(infs.fill)[1:33]), y = infs.fill[1,c(1:33)], type = 'l',
     ylim = c(-1,1), xlab = 'Harvest cycle', ylab = 'prcc')
for(i in 2:nrow(infs.fill)){
  lines(x = as.numeric(colnames(infs.fill)[1:33]), y = infs.fill[i,c(1:33)], col = i)
}

plot(x = as.numeric(colnames(ws.fill)[1:33]), y = ws.fill[1,c(1:33)], type = 'l',
     ylim = c(-1,1), xlab = 'Harvest cycle', ylab = 'prcc')
for(i in 2:nrow(ws.fill)){
  lines(x = as.numeric(colnames(ws.fill)[1:33]), y = ws.fill[i,c(1:33)], col = i)
}

plot(x = as.numeric(colnames(prev.fill)[1:33]), y = prev.fill[1,c(1:33)], type = 'l',
     ylim = c(-1,1), xlab = 'Harvest cycle', ylab = 'prcc')
for(i in 2:nrow(prev.fill)){
  lines(x = as.numeric(colnames(prev.fill)[1:33]), y = prev.fill[i,c(1:33)], col = i)
}
