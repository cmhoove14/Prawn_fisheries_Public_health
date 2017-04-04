#Assess sensitivity of infection outcomes, holding prawn dynamics constant to start ############
  source('snail_prawn_full.R') #Reset
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
