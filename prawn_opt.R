## Modeling prawn growth and mortality to maximize biomass

## To-do list:
  # Fit density-dependent parameters to data from Ranjeet & Kurup 2010 (or other Macrobrachium fisheries data)
    # Note: may also want to use net production data for fitting?
  # Investigate price-to-cost ratios in other locations for sensitivity analysis on maximum profit estimation

require(deSolve)

prawn_biomass = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
       
    Bp = (a/10)*(L/10)^b # Mean prawn mass, converting from length (mm) to weight (g)

    Bm = P*Bp # Total prawn biomass

    dLdt = k/(1+gam*Bm)*(linf - L) # Mean prawn length, using von Bertalanffy growth limited by density-dependent parameter gamma
    
    dPdt = -P*(mu*(Bp)^d + phi*Bm) # Prawn abundance, subject to size- and density-dependent mortality

    return(list(c(dPdt,dLdt)))
  }) 
} 

# Set initial values and parameters
# Stocking conditions: L = 25mm ~ W = 0.2g
nstart = c(P = 5000, L = 25)
time = seq(0, 365*2, 1)

parameters=c(
  a = 0.096868,        # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  b = 3.2944,          # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  gam = 1e-5,          # Density-dependent growth parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  mu = 0.00610958904,  # Prawn mortality at unit weight, from Lorenzen 1996 (pond aquaculture)
  d = -0.382,          # Exponential relationship of weight with mortality, from Lorenzen 1996 (pond aquaculture)
  phi = 5e-8,          # Density-dependent mortality parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  k = 0.00339726,      # Growth rate (mm/day), from Nwosu & Wolfi 2006 (M. vollenhovenii); alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333
  linf = 206           # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii)
)

# Economic parameters (price estimates from Tamil Nadu Agricultural University, http://agritech.tnau.ac.in/fishery/fish_freshwaterprawn.html)
p = 140                                           # Weighted average market price of prawns, in rupees/kg 
c = 600                                           # Cost of post-larvae, in rupees/1000 PL
delta = -log(1-0.1)/365                           # Discount rate, equivalent to 10%/year


# Run & plot
output = as.data.frame(ode(nstart,time,prawn_biomass,parameters))
output$B = ((parameters['a']*(output$L/10)^parameters['b'])/10)                    # Mean prawn biomass, transformed from length
output$Bt = output$B*output$P                                                      # Total prawn biomass
output$profit = p*(output$Bt/1000)*exp(-delta*(output$t)) - c*(nstart["P"]/1000)   # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 

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
x = c(1:50)
harvest = numeric(length(x))
harvest.m = numeric(length(x))
harvest.t = numeric(length(x))
for (i in x) {
  nstart.d = c(P = 1000*i, L = 25)
  output.d = as.data.frame(ode(nstart.d,time,prawn_biomass,parameters))
  output.d$B = ((parameters['a']*(output.d$L/10)^parameters['b'])/10)
  output.d$Bt = output.d$B*output.d$P
  harvest[i] = max(output.d$Bt)/1000
  harvest.m[i] = max(output.d$Bt)/output.d$Bt[output.d$time==1]
  harvest.t[i] = output.d$time[output.d$Bt == max(output.d$Bt)]
}
plot(x, harvest, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Harvest (kg/ha)',
     type = 'l', lwd = 2, main = 'Harvest by initial stocking density')
plot(x, harvest.m, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Biomass multiplier', 
     type = 'l', lwd = 2, main = 'Biomass multiplier by stocking density')
plot(x, harvest.t, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Harvest time (days)',
     type = 'l', lwd = 2, main = 'Harvest time by stocking density')

# Profit optimization
profit = p*harvest*exp(-delta*(harvest.t)) - c*x
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     type = 'l', lwd = 2, main = 'Estimated profit by stocking density')
  abline(v = which.max(profit), lty = 2, lwd = 2)


# # Fit density-dependent parameters to data from Ranjeet & Kurup 2010
# beta.params = function(mu, var) {
#   alpha = ((1 - mu)/var - 1/mu)*mu^2
#   beta = alpha*(1/mu - 1)
#   return(list(alpha = alpha, beta = beta))
# }
# 
# prawn.mle = function(params) {
#   parameters['gam'] = params[1]
#   parameters['phi'] = params[2]
#   
#   time = seq(0,30*8,1)
#   
#   nstart1 = c(P = 5000, L = 25)
#   output1 = as.data.frame(ode(nstart1, time, prawn_biomass, parameters))
#   mean.size1 = ((parameters['a']*(output1$L[output1$time == max(time)]/10)^parameters['b'])/10)
#   surv1 = output1$P[output1$time == max(time)] / output1$P[output1$time == 0]
#     
#   nstart2 = c(P = 10000, L = 25)
#   output2 = as.data.frame(ode(nstart2, time, prawn_biomass, parameters))
#   mean.size2 = ((parameters['a']*(output2$L[output2$time == max(time)]/10)^parameters['b'])/10)
#   surv2 = output2$P[output2$time == max(time)] / output2$P[output2$time == 0]
#     
#   nstart3 = c(P = 15000, L = 25)
#   output3 = as.data.frame(ode(nstart3, time, prawn_biomass, parameters))
#   mean.size3 = ((parameters['a']*(output3$L[output3$time == max(time)]/10)^parameters['b'])/10)
#   surv3 = output3$P[output3$time == max(time)] / output3$P[output3$time == 0]
#     
#   nstart4 = c(P = 25000, L = 25)
#   output4 = as.data.frame(ode(nstart4, time, prawn_biomass, parameters))
#   mean.size4 = ((parameters['a']*(output4$L[output4$time == max(time)]/10)^parameters['b'])/10)
#   surv4 = output4$P[output4$time == max(time)] / output4$P[output4$time == 0]
#   
#   lh.sz1 = dnorm(mean.size1, mean = 101.650, sd = 7.5067)
#   lh.sz2 = dnorm(mean.size2, mean = 87.730, sd = 8.1917)
#   lh.sz3 = dnorm(mean.size3, mean = 69.113, sd = 8.6865)
#   lh.sz4 = dnorm(mean.size4, mean = 55.483, sd = 3.9793)
#   
#   mu1 = 0.69440; var1 = 0.066525^2
#   lh.sv1 = dbeta(surv1, shape1 = beta.params(mu1, var1)$alpha, shape2 = beta.params(mu1, var1)$beta)
#   mu2 = 0.54533; var2 = 0.038724^2
#   lh.sv2 = dbeta(surv2, shape1 = beta.params(mu2, var2)$alpha, shape2 = beta.params(mu2, var2)$beta)
#   mu3 = 0.4269; var3 = 0.063707^2
#   lh.sv3 = dbeta(surv3, shape1 = beta.params(mu3, var3)$alpha, shape2 = beta.params(mu3, var3)$beta)
#   mu4 = 0.28210; var4 = 0.021198^2
#   lh.sv4 = dbeta(surv4, shape1 = beta.params(mu4, var4)$alpha, shape2 = beta.params(mu4, var4)$beta)
#   
#   negLL = -(log(lh.sz1)+log(lh.sz2)+log(lh.sz3)+log(lh.sz4)+log(lh.sv1)+log(lh.sv2)+log(lh.sv3)+log(lh.sv4))
#   return(negLL)
# }
# 
# prawn.mle(c(1e-6, 4e8))
# 
# ## DOESN'T WORK - ODE solver runs into problems during optimization
# density.params.fit = optim(par=c(1e-6, 4e8), prawn.mle, method = 'Nelder-Mead')


