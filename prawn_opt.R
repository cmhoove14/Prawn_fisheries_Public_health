## Modeling prawn growth and mortality to maximize biomass

## To-do list:
  # Fit density-dependent parameters to data from Ranjeet & Kurup 2010 (or other Macrobrachium fisheries data)
    # Note: update distributions used for fitting - survival shouldn't be normal (use beta distribution?)
  # Look into other parameter estimates: M. rosenbergii vs. M. vollenhovenii, extensive vs. intensive aquaculture
  # Try to optimize profit based on stocking density of post-larvae
    # Objective f'n: profit = p*H(SD)*e^(-d*T(SD)) - c*SD

require(deSolve)

prawn_biomass = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
       
    Bp = (a*(L/10)^b)/10 # Mean prawn mass, converting from length (mm) to weight (g)

    Bm = P*Bp # Total prawn biomass

    dLdt = k/(1+gam*Bm)*(linf - L) # Mean prawn length, using von Bertalanffy growth limited by density-dependent parameter gamma
    
    dPdt = -P*(mu*L^d + Bm/phi) # Prawn abundance, subject to size- and density-dependent mortality

    return(list(c(dPdt,dLdt)))
  }) 
} 

# Set initial values and parameters
nstart = c(P = 15000, L = 25)
time = seq(0, 365*2, 1)

parameters=c(
  a = 0.096868,     # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (growout phase)
  b = 3.2944,       # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (growout phase)
  gam = 1e-6,       # Density-dependent growth parameter (based on biomass per hectare); NEEDS TO BE FIT
  mu = 0.006136986, # Natural prawn mortality rate, from Nwosu & Wolfi 2006
  d = -0.25,        # Exponential relationship of size with mortality; no source
  phi = 40000000,   # Density-dependent mortality parameter (based on biomass per hectare); NEEDS TO BE FIT
  k = 0.00339726,   # Growth rate (mm/day), from Nwosu & Wolfi 2006
  linf = 206        # Max length (mm), from Nwosu & Wolfi 2006
)


# Run & plot
output = as.data.frame(ode(nstart,time,prawn_biomass,parameters))
output$B = ((parameters['a']*(output$L/10)^parameters['b'])/10)    # Mean prawn biomass, transformed from length
output$Bt = output$B*output$P                                      # Total prawn biomass

start.mass.kg = output$Bt[output$time==1]/1000
harvest.mass.kg = max(output$Bt)/1000
harvest.size = output$B[output$Bt==max(output$Bt)]  
harvest.time = output$time[output$Bt==max(output$Bt)]

plot(x = output$time, y = output$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(output$time)), ylim = c(0, max(output$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart[2]), ' mm)', sep = ''))
  lines(x = output$time, y = output$Bt/1000, col = 'blue', lwd=2)
  lines(x = output$time, y = output$B, col = 'purple', lwd=2, lty=2)
  lines(x = output$time, y = output$L, col = 'green', lwd=2)
  abline(v = harvest.time, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5)
  legend('bottomright', legend = c(paste('Starting mass =', round(start.mass.kg), 'kg', sep = ' '),
                           paste('Total harvest mass =', round(harvest.mass.kg), 'kg', sep = ' '), 
                           paste('Mean harvest mass =', round(harvest.size), 'g', sep = ' '), 
                           paste('Time of harvest =', round(harvest.time), 'days', sep = ' ')), cex=0.5)
  
# Plot amount of harvest, biomass multiplier, optimal harvest time against initial stocking density
x = c(0:50)
harvest = numeric(length(x))
harvest.m = numeric(length(x))
harvest.t = numeric(length(x))
for (i in x) {
  nstart.d = c(P = 1000*i, L = 25)
  output.d = as.data.frame(ode(nstart.d,time,prawn_biomass,parameters))
  output.d$B = ((parameters['a']*(output.d$L/10)^parameters['b'])/10)
  output.d$Bt = output.d$B*output.d$P
  harvest[i+1] = max(output.d$Bt)/1000
  harvest.m[i+1] = max(output.d$Bt)/output.d$Bt[output.d$time==1]
  harvest.t[i+1] = output.d$time[output.d$Bt == max(output.d$Bt)]
}
plot(x, harvest, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Harvest (kg/ha)',
     type = 'l', lwd = 2, main = 'Harvest by initial stocking density')
plot(x[-1], harvest.m[-1], col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Biomass multiplier', 
     type = 'l', lwd = 2, main = 'Biomass multiplier by stocking density')
plot(x[-1], harvest.t[-1], col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Harvest time (days)',
     type = 'l', lwd = 2, main = 'Harvest time by stocking density')


# Fit density-dependent parameters to data from Ranjeet & Kurup 2010
prawn.mle = function(gam, phi) {
  parameters['gam'] = gam
  parameters['phi'] = phi
  
  time = seq(0,30*8,1)
  
  nstart1 = c(P = 5000, L = 25)
  output1 = as.data.frame(ode(nstart1,time,prawn_biomass,parameters))
  mean.size1 = ((parameters['a']*(output1$L[output1$time == max(time)]/10)^parameters['b'])/10)
  surv1 = output1$P[output1$time == 0] / output1$P[output1$time == max(time)]
    
  nstart2 = c(P = 10000, L = 25)
  output2=as.data.frame(ode(nstart2,time,prawn_biomass,parameters))
  mean.size2 = ((parameters['a']*(output2$L[output2$time == max(time)]/10)^parameters['b'])/10)
  surv2 = output2$P[output2$time == 0] / output2$P[output2$time == max(time)]
    
  nstart3 = c(P = 15000, L = 25)
  output3=as.data.frame(ode(nstart3,time,prawn_biomass,parameters))
  mean.size3 = ((parameters['a']*(output3$L[output3$time == max(time)]/10)^parameters['b'])/10)
  surv3 = output3$P[output3$time == 0] / output3$P[output3$time == max(time)]
    
  nstart4 = c(P = 25000, L = 25)
  output4=as.data.frame(ode(nstart4,time,prawn_biomass,parameters))
  mean.size4 = ((parameters['a']*(output4$L[output4$time == max(time)]/10)^parameters['b'])/10)
  surv4 = output4$P[output4$time == 0] / output4$P[output4$time == max(time)]
  
  ll1 = dnorm(mean.size1, mean = 101.650, sd = 7.5067)
  ll2 = dnorm(surv1, mean = 0.69440, sd = 0.066525)
  ll3 = dnorm(mean.size2, mean = 87.730, sd = 8.1917)
  ll4 = dnorm(surv2, mean = 0.54533, sd = 0.038724)
  ll5 = dnorm(mean.size3, mean = 69.113, sd = 8.6865)
  ll6 = dnorm(surv3, mean = 0.4269, sd = 0.063707)
  ll7 = dnorm(mean.size4, mean = 55.483, sd = 3.9793)
  ll8 = dnorm(surv4, mean = 0.28210, sd = 0.021198)

  LL.fin = ll1*ll2*ll3*ll4*ll5*ll6*ll7*ll8
  negLL = -log(LL.fin+1)
  return(negLL)
}

prawn.optim(params = c(1e-7, 3e8))

op.df = data.frame('gams' = seq(1e-7, 1e-5, 1e-7),
                   'phis' = seq(1e7, 1e9, 1e7),
                   'negLL' = 0) 
for(i in 1:nrow(op.df)){
  params = c(op.df[i,1], op.df[i,2])
  op.df[i,3] = prawn.optim(params)
}

op.bm<-optim(par=c(1e-6, 4e8), prawn.optim, method = 'Nelder-Mead', control = list(fnscale = -1))

