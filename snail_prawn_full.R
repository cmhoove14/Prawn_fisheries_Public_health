#### Full snail-prawn model including epidemiological, predation, and aquaculture components

## To-do list:
  # Investigate the long-term stability of snail elimination, particularly at low attack rates
  # Incorporate seasonality into model
  # Sensitivity analysis, particularly for prawn attack rate

require(deSolve)

snail_prawn_model = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S1=n[1]
    S2=n[2]
    S3=n[3]
    E1=n[4]
    E2=n[5]
    E3=n[6]
    I2=n[7]
    I3=n[8]
    W=n[9]
    P=n[10]
    L=n[11]
    
    N1 = S1+E1     # Total snails of size class 1
    N2 = S2+E2+I2  # Total snails of size class 2
    N3 = S3+E3+I3  # Total snails of size class 3
    N = N1+N2+N3   # Total number of snails
    
    # Miracidial density per person as a function of mean worm burden (W) and miracidial shedding rate (m)
    M = m*0.5*W
    
    # Mean and total prawn biomass, converting from length (mm) to weight (g)
    Bm.p = (a.p*(L/10)^b.p)/10  
    Bm.t = P*Bm.p
    
    # Mean snail biomass in each size class, converting from length (cm) to weight (g)
    Bm.n1 = a.s*0.4^b.s  # 4mm size class
    Bm.n2 = a.s*0.8^b.s  # 8mm size class
    Bm.n3 = a.s*1.2^b.s  # 12mm size class
    
    # Prawn-to-snail mass ratios for each size class
    Bm.r1 = Bm.p / Bm.n1
    Bm.r2 = Bm.p / Bm.n2
    Bm.r3 = Bm.p / Bm.n3
    
    # Attack rates for each size class as a function of biomass ratio
    alpha1 = ar*Bm.r1
    alpha2 = ar*Bm.r2
    alpha3 = ar*Bm.r3
    
    # Adjusted attack rates, accounting for area of interest and limiting factors in the wild
    alpha_star1 = alpha1*s/A
    alpha_star2 = alpha2*s/A
    alpha_star3 = alpha3*s/A
    
    # Handling times for each size class as a function of biomass ratio
    handle1 = 1/(th*Bm.r1)
    handle2 = 1/(th*Bm.r2)
    handle3 = 1/(th*Bm.r3)
    
    # Functional responses for each size class
    psi1 = (P*alpha_star1) / (1+sum(alpha_star1*handle1*N1, alpha_star2*handle2*N2, alpha_star3*handle3*N3))
    psi2 = (P*alpha_star2) / (1+sum(alpha_star1*handle1*N1, alpha_star2*handle2*N2, alpha_star3*handle3*N3))
    psi3 = (P*alpha_star3) / (1+sum(alpha_star1*handle1*N1, alpha_star2*handle2*N2, alpha_star3*handle3*N3))
    
    ## Model equations:
    
    # Susceptible snails, class 1
    dS1dt = f*(1-N/(Kn*A))*(S2+S3+z*(E2+E3)) - muN1*S1 - psi1*S1 - g1*S1 - beta*M*H*S1
    
    # Susceptible snails, class 2
    dS2dt = g1*S1 - muN2*S2 - psi2*S2 - g2*S2 - beta*M*H*S2
    
    # Susceptible snails, class 3
    dS3dt = g2*S2 - muN3*S3 - psi3*S3 - beta*M*H*S3
    
    # Exposed (prepatent) snails, class 1
    dE1dt = beta*M*H*S1 - muN1*E1 - psi1*E1 - g1*E1 - sigma*X*E1
    
    # Exposed (prepatent) snails, class 2
    dE2dt = beta*M*H*S2 + g1*E1 - muN2*E2 - psi2*E2 - g2*E2 - sigma*E2
    
    # Exposed (prepatent) snails, class 3
    dE3dt = beta*M*H*S3 + g2*E2 - muN3*E3 - psi3*E3 - sigma*E3
    
    # Infectious (shedding) snails, class 2
    dI2dt = sigma*X*E1 + sigma*(1-rho)*E2 - (muN2 + muI)*I2 - psi2*I2 - g2*I2
    
    # Infectious (shedding) snails, class 3
    dI3dt = sigma*rho*E2 + sigma*E3 + g2*I2 - (muN3 + muI)*I3 - psi3*I3
    
    # Mean human worm burden
    dWdt = lambda*I2/A + theta*lambda*I3/A - (muW + muH)*W
    
    # Prawn abundance
    dPdt = -P*(muP*L^d + Bm.t/phi)
    
    # Mean prawn length (mm)
    dLdt = k/(1+gam*Bm.t)*(linf - L)
    
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI2dt, dI3dt, dWdt, dPdt, dLdt)))
  }) 
} 


## Set initial values and parameters
#  Default settings: area = 10000, P = 11000/ha, L = 25 (0.2g post-larvae)
#  Gates settings: area = 20000; P = 500, 1000, 2500, 5000, 10000/ha; L = 67 or 100 (5g or 20g adults)
area = 20000
nstart = c(S1 = 0.144*area, S2 = 0.002*area, S3 = 0, E1 = 6.57*area, E2 = 2.61*area, E3 = 0.843*area, 
           I2 = 0.640*area, I3 = 0.329*area, W = 74, P = 500*area/10000, L = 100)
time = seq(0, 365*2, 1)

parameters=c(
  # Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1000,          # Human population at site of interest
  
  # Snail reproductive parameters
  f = 0.26,          # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection - more like a recruitment rate); 
                     #   adjusted from Sokolow et al. 2015 to account for the fact that class 1 snails are assumed to be juveniles and therefore don't reproduce
  Kn = 50,           # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.5,           # Fraction of prepatent snails that can reproduce, from Sokolow et al. 2015
  
  # Snail growth parameters
  a.s = 0.187178454, # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
  b.s = 2.536764792, # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
  g1 = 1/37,         # Growth rate of small snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  g2 = 1/62,         # Growth rate of medium snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  
  # Snail mortality parameters
  muN1 = 1/40,       # Natural mortality rate of small snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN2 = 1/50,       # Natural mortality rate of medium snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN3 = 1/60,       # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 50 days)
  muI = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  # Prawn growth parameters
  a.p = 0.096868,    # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  b.p = 3.2944,      # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  k = 0.00339726,    # Growth coefficient, from Nwosu & Wolfi 2006 (M. vollenhovenii); alternate value for M. rosenbergii, from Sampaio & Wagner 1996: 0.0104333333
  linf = 206,        # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii)
  gam = 1e-6,        # Density-dependent growth parameter (based on biomass per hectare); not yet fit
  
  # Prawn mortality parameters
  muP = 0.006136986, # Natural prawn mortality rate, from Nwosu & Wolfi 2006 (M. vollenhovenii)
  d = -0.25,         # Exponential parameter relating size with mortality; no source
  phi = 2e7,         # Density-dependent mortality parameter (based on biomass per hectare); not yet fit
  
  # Predation parameters
  ar = 0.037192,     # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
  th = 0.40450,      # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  s = 1/50,          # Scale factor limiting prawn attack rate in the wild (vs. lab conditions);
                     #   decreased from 1/10 to 1/50 in size class model to match expected elimination threshold
  
  # Infection parameters
  beta = 8e-5,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day); adjusted from Sokolow et al. 2015 (original value: 4e-6)
  m = 0.8,           # Miracidial shedding rate per adult female worm divided by miracidial mortality; from Sokolow et al. 2015
  sigma = 1/30,      # Latent period for exposed snails (infectious snails/exposed snail/day); adjusted from Sokolow et al. 2015 (original value: 1/50)
  X = 0,             # Fraction of exposed small snails that convert directly to medium shedding (ignored for now)
  rho = 0,           # Fraction of exposed medium snails that convert directly to large shedding (ignored for now)
  lambda = 0.05,     # Snail-to-human infection probability scaled to 1 m^2 (composite including cercarial shedding, mortality, infection, survival to patency); 
                     #   adjusted from Sokolow et al. 2015 (original value: 0.1)
  theta = 2,         # Scale factor describing increase in cercarial shedding rate in larger snails; from Chu & Dawood 1970 (estimated to be between 2 and 10)
  
  # Schisto mortality parameters
  muW = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  muH = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)



## Run model and assess outcomes over one or more aquaculture cycles

# Start with single run to calculate the time to optimal harvest and other fisheries statistics
output = as.data.frame(ode(nstart, time, snail_prawn_model, parameters))
output$B = ((parameters['a.p']*(output$L/10)^parameters['b.p'])/10)       # Mean prawn biomass, transformed from length using allometric equation
output$Bt = output$B*output$P                                             # Total prawn biomass

start.mass.kg = output$Bt[output$time == 0]/1000                          # Total prawn biomass at beginning of cycle, in kg
harvest.mass.kg = max(output$Bt)/1000                                     # Total prawn biomass at end of cycle, in kg
harvest.size = output$B[output$Bt == max(output$Bt)]                      # Average prawn biomass at harvest, in g
harvest.time = output$time[output$Bt == max(output$Bt)]                   # Time to optimal harvest, in days

# Plot prawn aquaculture dynamics
plot(x = output$time, y = output$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(output$time)), ylim = c(0, max(output$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart[11]), ' mm)', sep = ''))
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

# Set the desired number of aquaculture cycles, and harvest time if not using optimum
# Default settings: harvest.time = optimal, ncycles = 15 (~10 years @ 8 months/cycle)
# Gates settings: harvest.time = 4 months, ncycles = 3 (1 year)
harvest.time = 365/3
ncycles = 3
pzq.delay = 0
nstart.lt = nstart
nstart.lt['P'] = 0
time.lt = seq(0, pzq.delay + harvest.time*ncycles, 1)

# Define stocking events at the beginning of each cycle
stocking = data.frame(var = rep(c('P', 'L'), ncycles),
                      time = rep(seq(pzq.delay, pzq.delay + harvest.time*(ncycles-1), harvest.time), each = 2),
                      value = rep(c(nstart['P'], nstart['L']), ncycles),
                      method = rep('rep', ncycles*2))

# Define PZQ administration event
# Note: if administering before start of stocking, set time = 0 and instead set pzq.delay above
#   to time between MDA and first stocking
pzq = data.frame(var = 'W', time = 30, value = 2, method = 'rep')

# Run model and calculate outcomes of interest
output.lt = as.data.frame(ode(nstart.lt, time.lt, snail_prawn_model, parameters,
                              events = list(data = rbind(stocking, pzq))))

output.lt$S.t = output.lt$S1 + output.lt$S2 + output.lt$S3                      # Total susceptible snails
output.lt$E.t = output.lt$E1 + output.lt$E2 + output.lt$E3                      # Total exposed snails 
output.lt$I.t = output.lt$I2 + output.lt$I3                                     # Total infected snails
output.lt$N1.t = output.lt$S1 + output.lt$E1                                    # Total snails of size class 1
output.lt$N2.t = output.lt$S2 + output.lt$E2 + output.lt$I2                     # Total snails of size class 2
output.lt$N3.t = output.lt$S3 + output.lt$E3 + output.lt$I3                     # Total snails of size class 3
output.lt$N.t = output.lt$S.t + output.lt$E.t + output.lt$I.t                   # Total snails
output.lt$prev = pnbinom(2, size = 0.25, mu = output.lt$W, lower.tail = FALSE)  # Estimated prevalence, using a negative binomial dist. with k = 0.25
output.lt$B = ((parameters['a.p']*(output.lt$L/10)^parameters['b.p'])/10)       # Mean prawn biomass, transformed from length using allometric equation
output.lt$Bt = output.lt$B*output.lt$P                                          # Total prawn biomass

# Plot snail dynamics by size class
plot(x = output.lt$time, y = output.lt$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output.lt$N.t)),
     main = 'Snail Size Classes')
lines(output.lt$time, output.lt$N1.t, col = 'green', lwd = 2)
lines(output.lt$time, output.lt$N2.t, col = 'blue', lwd = 2)
lines(output.lt$time, output.lt$N3.t, col = 'red', lwd = 2)
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



## Assess single-cycle outcomes over multiple stocking densities
#  NOTE: set parameters and initial values first!

# Run the model over the desired range of stocking densities (in thousands of post-larvae)
x = c(0:30)
nstart.sd = nstart
time.sd = seq(0, 500, 1)
harvest = numeric(length(x))
harvest.t = numeric(length(x))
snails = numeric(length(x))
worms = numeric(length(x))
for (i in x) {
  nstart.sd['P'] = 1000*i
  output.sd = as.data.frame(ode(nstart.sd, time.sd, snail_prawn_model, parameters,
                                events = list(data = pzq)))
  output.sd$S.t = output.sd$S1 + output.sd$S2 + output.sd$S3
  output.sd$E.t = output.sd$E1 + output.sd$E2 + output.sd$E3
  output.sd$I.t = output.sd$I2 + output.sd$I3
  output.sd$Inf.t = output.sd$E.t + output.sd$I.t
  output.sd$N.t = output.sd$S.t + output.sd$E.t + output.sd$I.t
  output.sd$B = ((parameters['a.p']*(output.sd$L/10)^parameters['b.p'])/10)
  output.sd$Bt = output.sd$B*output.sd$P
  harvest[i+1] = max(output.sd$Bt)/1000
  ht.tmp = output.sd$time[output.sd$Bt == max(output.sd$Bt)]
  harvest.t[i+1] = ifelse(ht.tmp != 0, ht.tmp, max(output.sd$time))
  snails[i+1] = output.sd$N.t[output.sd$time == harvest.t[i+1]]
  worms[i+1] = output.sd$W[output.sd$time == harvest.t[i+1]]
}

# Price estimates for profit calculation from Tamil Nadu Agricultural University, http://agritech.tnau.ac.in/fishery/fish_freshwaterprawn.html
p = 140                                           # Weighted average market price of prawns, in rupees/kg 
c = 600                                           # Cost of post-larvae, in rupees/1000 PL
delta = -log(1-0.1)/365                           # Discount rate, equivalent to 10%/year
profit = p*harvest*exp(-delta*(harvest.t)) - c*x  # Profit function, in terms of revenue (discounted by time to harvest) minus stocking costs 

# Plot estimated profit and time to harvest over stocking density
# (after one aquaculture cycle)
par(mar = c(5,5,6,5))
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     type = 'l', lwd = 2, main = 'Profit and time to harvest \n by stocking density')
abline(v = which.max(profit)-1, lty = 2, lwd = 2)
par(new = T)
plot(x, harvest.t, col = 'blue', axes = F, xlab = NA, ylab = NA, ylim = c(0, max(harvest.t[-1])), type = 'l', lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Harvest time (days)')
legend('topright', legend = c('Profit', 'Harvest time'), lty = 1, col = c('green', 'blue'), cex = 0.7)

# Plot estimated profit and snail abundance over stocking density
# (after one aquaculture cycle)
par(mar = c(5,5,6,5))
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     type = 'l', lwd = 2, main = 'Profit and snail abundance \n by stocking density')
abline(v = which.max(profit)-1, lty = 2, lwd = 2)
par(new = T)
plot(x, snails, col = 'orange', axes = F, xlab = NA, ylab = NA, ylim = c(0, max(snails)), type = 'l', lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Number of snails')
legend('topright', legend = c('Profit', 'Snails'), lty = 1, col = c('green', 'orange'), cex = 0.7)

# Plot estimated profit and mean worm burden over stocking density
# (after one aquaculture cycle)
par(mar = c(5,5,6,5))
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     type = 'l', lwd = 2, main = 'Profit and worm burden (after PZQ) \n by stocking density')
abline(v = which.max(profit)-1, lty = 2, lwd = 2)
par(new = T)
plot(x, worms, col = 'red', axes = F, xlab = NA, ylab = NA, ylim = c(0, max(worms)), type = 'l', lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Mean worm burden')
legend('topright', legend = c('Profit', 'Mean worm burden'), lty = 1, col = c('green', 'red'), cex = 0.7)




