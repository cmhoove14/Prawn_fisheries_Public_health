## Modeling prawn growth and mortality to maximize biomass

## To-do list: #########
  # Make cost of initial prawn stocking a function of their size (unless we want to fix initial stocking size)
  # Make price of harvested prawns a function of their size such that larger prawns are more valuable
  # Fit density-dependent parameters to data from Ranjeet & Kurup 2010 (or other Macrobrachium fisheries data)
    # Note: may also want to use net production data for fitting?
  # Investigate price-to-cost ratios in other locations for sensitivity analysis on maximum profit estimation

require(plotly)

source('Prawn_aquaculture/prawn_aquaculture_mod.R')

# Run & plot typical cycle ######
output = as.data.frame(ode(nstart.p,t.p,prawn_biomass,par.aqua))
output$B = ((par.aqua['a.p']*(output$L/10)^par.aqua['b.p'])/10)                    # Mean prawn biomass, transformed from length
output$Bt = output$B*output$P                                                      # Total prawn biomass
output$profit = p*(output$Bt/1000)*exp(-delta*(output$t)) - c*(nstart.p["P"]/1000)   # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 

start.mass.kg = output$Bt[output$time==0]/1000
harvest.mass.kg = max(output$Bt)/1000
harvest.size = output$B[output$Bt==max(output$Bt)]  
harvest.P = output$P[output$Bt==max(output$Bt)]  
harvest.time = output$time[output$Bt==max(output$Bt)]

plot(x = output$time, y = output$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(output$time)), ylim = c(0, max(output$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.p[2]), ' mm)', sep = ''))
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
  
# Plot amount of harvest, biomass multiplier, optimal harvest time against initial stocking density #######
x = c(1:50)
harvest = numeric(length(x))
harvest.m = numeric(length(x))
harvest.t = numeric(length(x))
for (i in x) {
  nstart.d = c(P = 1000*i, L = 25)
  output.d = as.data.frame(ode(nstart.d,t.p,prawn_biomass,par.aqua))
  output.d$B = ((par.aqua['a.p']*(output.d$L/10)^par.aqua['b.p'])/10)
  output.d$Bt = output.d$B*output.d$P
  harvest[i] = max(output.d$Bt)/1000
  harvest.m[i] = max(output.d$Bt)/output.d$Bt[output.d$time==0]
  harvest.t[i] = output.d$time[output.d$Bt == max(output.d$Bt)]
}
plot(x, harvest, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Harvest (kg/ha)',
     type = 'l', lwd = 2, main = 'Harvest by initial stocking density')
plot(x, harvest.m, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Biomass multiplier', 
     type = 'l', lwd = 2, main = 'Biomass multiplier by stocking density')
plot(x, harvest.t, col = 'blue', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Harvest time (days)',
     type = 'l', lwd = 2, main = 'Harvest time by stocking density')

# Profit optimization #######
profit = p*harvest*exp(-delta*(harvest.t)) - c*x
ROI = (p*harvest*exp(-delta*(harvest.t)) - c*x) / (c*x)

plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     type = 'l', lwd = 2, main = 'Estimated profit by stocking density')
  abline(v = which.max(profit), lty = 2, lwd = 2)
  
plot(x, ROI, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Return on investment',
     type = 'l', lwd = 2, main = 'Estimated ROI by stocking density')
  abline(v = which.max(ROI), lty = 2, lwd = 2)

# Check optimization wrt to starting prawn length ###########
l = c(1:50)
  harvest.l = numeric(length(l))
  harvest.l.m = numeric(length(l))
  harvest.l.t = numeric(length(l))
  for (i in l) {
    nstart.l = c(P = 5000, L = i)
    output.l = as.data.frame(ode(nstart.l,t.p,prawn_biomass,par.aqua))
    output.l$B = ((par.aqua['a.p']*(output.l$L/10)^par.aqua['b.p'])/10)
    output.l$Bt = output.l$B*output.l$P
    harvest.l[i] = max(output.l$Bt)/1000
    harvest.l.m[i] = max(output.l$Bt)/output.l$Bt[output.l$time==0]
    harvest.l.t[i] = output.l$time[output.l$Bt == max(output.l$Bt)]
  }
  plot(l, harvest.l, col = 'blue', xlab = 'Prawn length (mm)', ylab = 'Harvest (kg/ha)',
       type = 'l', lwd = 2, main = 'Harvest by initial prawn length')
  plot(l, harvest.l.m, col = 'blue', xlab = 'Prawn length (mm)', ylab = 'Biomass multiplier', 
       type = 'l', lwd = 2, main = 'Biomass multiplier by initial prawn length')
  plot(l, harvest.l.t, col = 'blue', xlab = 'Prawn length (mm)', ylab = 'Harvest time (days)',
       type = 'l', lwd = 2, main = 'Harvest time by initial prawn length')
  
profit.l = p*harvest.l*exp(-delta*(harvest.l.t)) - c*x
  plot(l, profit.l, col = 'green', xlab = 'Prawn length (mm)', ylab = 'Profit (rupees/ha)',
       type = 'l', lwd = 2, main = 'Estimated profit by initial prawn length')
    abline(v = which.max(profit.l), lty = 2, lwd = 2)
    
ROI.l = (p*harvest.l*exp(-delta*(harvest.l.t)) - c*x) / (c*x)
    
  plot(x, ROI, col = 'green', xlab = 'Prawn length (mm)', ylab = 'Return on investment',
       type = 'l', lwd = 2, main = 'Estimated ROI by initial prawn length')
    abline(v = which.max(ROI), lty = 2, lwd = 2)
    

#plot biomass trajectory with different starting lengths and densities ###########
  bm.sim = expand.grid(l0 = c(5,15,25), p0 = c(1000, 5000, 10000)) 
  
  bm.arr = array(data = NA, dim = c(length(t.p), ncol(output), nrow(bm.sim)))
  
  for(k in 1:nrow(bm.sim)){
    start = c(P = bm.sim[k,2], L = bm.sim[k,1])
    op = as.data.frame(ode(start,t.p,prawn_biomass,par.aqua))
    op$B = ((par.aqua['a.p']*(op$L/10)^par.aqua['b.p'])/10)   # Mean prawn biomass, transformed from length
    op$Bt = op$B*op$P/1000
    op$profit = p*(max(op$Bt)/1000)*exp(-delta*(op$time[op$Bt==max(op$Bt)])) - c*(start["P"]/1000)  
    
    bm.arr[ , , k] = as.matrix(op)
  }
  
  plot(bm.arr[ , 1, 1], bm.arr[, 5, 1], type = 'l', lwd = 2, lty = 3, col = 2,
       ylim = c(0,max(bm.arr[ , 5, ])), xlab = 'time', ylab = 'prawn biomass (kg)')
    lines(bm.arr[ , 1, 2], bm.arr[, 5, 2], lwd = 2, lty = 3, col = 3)
    lines(bm.arr[ , 1, 3], bm.arr[, 5, 3], lwd = 2, lty = 3, col = 4)
    lines(bm.arr[ , 1, 4], bm.arr[, 5, 4], lwd = 2, lty = 2, col = 2)
    lines(bm.arr[ , 1, 5], bm.arr[, 5, 5], lwd = 2, lty = 2, col = 3)
    lines(bm.arr[ , 1, 6], bm.arr[, 5, 6], lwd = 2, lty = 2, col = 4)
    lines(bm.arr[ , 1, 7], bm.arr[, 5, 7], lwd = 2, lty = 1, col = 2)
    lines(bm.arr[ , 1, 8], bm.arr[, 5, 8], lwd = 2, lty = 1, col = 3)
    lines(bm.arr[ , 1, 9], bm.arr[, 5, 9], lwd = 2, lty = 1, col = 4)
    
  
#Optimize by length and stocking density for M. volenhovenii simultaneously ######
opt.df = expand.grid(L_nought = c(1:50), P_nought = seq(500, 50000, 500)) 
  opt.df$h.t = 0
  opt.df$h.bm = 0
  opt.df$p.bm = 0
  opt.df$p.surv = 0
  opt.df$profit = 0
  opt.df$ROI = 0

  par.aqua['k'] = 0.00339726       # Growth rate (mm/day), from Nwosu & Wolfi 2006 (M. vollenhovenii)
  
  for(k in 1:nrow(opt.df)){
    start = c(P = opt.df[k,2], L = opt.df[k,1])
    op = as.data.frame(ode(start,t.p,prawn_biomass,par.aqua))
    op$B = ((par.aqua['a.p']*(op$L/10)^par.aqua['b.p'])/10)                    # Mean prawn biomass, transformed from length
    
    opt.df[k,3] = op$time[op$B*op$P==max(op$B*op$P)]  # harvest time
    opt.df[k,4] = max(op$B*op$P)                      # Total prawn biomass
    opt.df[k,5] = op$B[op$B*op$P==max(op$B*op$P)]     # Prawn size at harvest
    opt.df[k,6] = op$P[op$B*op$P==max(op$B*op$P)] / opt.df[k,2]   # Prawn % survival at harvest
    opt.df[k,7] = p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"]/1000)     #Profit
    opt.df[k,8] = (p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"]/1000)) / (c*(start["P"]/1000))#ROI  
    
    if(k %% 1000==0) print(k)
  }


#contour plots with no restrictions on harvest time, harvest size, or profit #########
plot_ly(opt.df, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')
plot_ly(opt.df, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')
plot_ly(opt.df, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')

#restrict to runs in which harvest time was <1 yearand profit was positive ##########
opt.df.pos = opt.df[which(opt.df$profit > 0 & opt.df$h.t <= 365),]

plot_ly(opt.df.pos, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')

plot_ly(opt.df.pos, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')
plot_ly(opt.df.pos, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')

#Restrict to runs in which initial stocking size is restricted ####################
#Bigger prawns at stocking inevitable leads to higher profits. 
#What's the max initial size feasible for stocking?
#say it's 25mm
opt.df.pos.25 = opt.df.pos[which(opt.df.pos$L_nought <= 25),]
plot_ly(opt.df.pos.25, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')

plot_ly(opt.df.pos.25, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')
plot_ly(opt.df.pos.25, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')

opt.df.pos.25$P_nought[opt.df.pos.25$profit == max(opt.df.pos.25$profit)]

#Restrict to runs in which sufficient marketable size is achieved ##############
#Also want to make sure that prawns are actually big enough to be marketable when harvested
#Restrict mean prawn biomass to be >= 40g when harvested
opt.df.pos.25.40g = opt.df.pos.25[which(opt.df.pos.25$p.bm >= 40),]

plot_ly(opt.df.pos.25.40g, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')

plot_ly(opt.df.pos.25.40g, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')
plot_ly(opt.df.pos.25.40g, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')

vol.opt = max(opt.df.pos.25.40g$ROI)

vol.l0 = opt.df.pos.25.40g$L_nought[opt.df.pos.25.40g$ROI == vol.opt]  # volenhovenii optimal starting length
vol.p0 = opt.df.pos.25.40g$P_nought[opt.df.pos.25.40g$ROI == vol.opt]  # volenhovenii optimal starging density
vol.ht = opt.df.pos.25.40g$h.t[opt.df.pos.25.40g$ROI == vol.opt]       # volenhovenii optimal harvest time
vol.hbm = opt.df.pos.25.40g$h.bm[opt.df.pos.25.40g$ROI == vol.opt]     # volenhovenii harvest biomass at optimal harvest
vol.pbm = opt.df.pos.25.40g$p.bm[opt.df.pos.25.40g$ROI == vol.opt]     # volenhovenii prawn mass at optimal harvest


#Optimize by length and stocking density for M. rosenbergii simultaneously ######
opt.dfr = expand.grid(L_nought = c(1:50), P_nought = seq(500, 50000, 500)) 
  opt.dfr$h.t = 0
  opt.dfr$h.bm = 0
  opt.dfr$p.bm = 0
  opt.dfr$p.surv = 0
  opt.dfr$profit = 0
  opt.dfr$ROI = 0

  par.aqua['k'] = 0.0104333333  # alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333

for(k in 1:nrow(opt.dfr)){
  start = c(P = opt.dfr[k,2], L = opt.dfr[k,1])
  op = as.data.frame(ode(start,t.p,prawn_biomass,par.aqua))
  op$B = ((par.aqua['a.p']*(op$L/10)^par.aqua['b.p'])/10)                    # Mean prawn biomass, transformed from length
  
  opt.dfr[k,3] = op$time[op$B*op$P==max(op$B*op$P)]  # harvest time
  opt.dfr[k,4] = max(op$B*op$P)                      # Total prawn biomass
  opt.dfr[k,5] = op$B[op$B*op$P==max(op$B*op$P)]     # Prawn size at harvest
  opt.dfr[k,6] = op$P[op$B*op$P==max(op$B*op$P)] / opt.dfr[k,2]   # Prawn % survival at harvest
  opt.dfr[k,7] = p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"]/1000)     #Profit
  opt.dfr[k,8] = (p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"]/1000)) / (c*(start["P"]/1000))#ROI  
  
  if(k %% 100==0) print(k)
}

#contour plots with no restrictions on harvest time, harvest size, or profit ############
  plot_ly(opt.dfr, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')
  plot_ly(opt.dfr, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')
  plot_ly(opt.dfr, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')

#restrict to runs in which harvest time was <1 year and profit was positive##########
opt.dfr.pos = opt.dfr[which(opt.dfr$profit > 0 & opt.dfr$h.t <= 365),]

  plot_ly(opt.dfr.pos, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')
  
  plot_ly(opt.dfr.pos, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')
  plot_ly(opt.dfr.pos, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')

#Restrict to runs in which initial stocking size is restricted ####################
#Bigger prawns at stocking inevitable leads to higher profits. 
#What's the max initial size feasible for stocking?
#say it's 25mm
opt.dfr.pos.25 = opt.dfr.pos[which(opt.dfr.pos$L_nought <= 25),]
  
  plot_ly(opt.dfr.pos.25, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')
  
  plot_ly(opt.dfr.pos.25, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')
  plot_ly(opt.dfr.pos.25, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')
  
#Restrict to runs in which sufficient marketable size is achieved ##############
#Also want to make sure that prawns are actually big enough to be marketable when harvested
#Restrict mean prawn biomass to be >= 40g when harvested
opt.dfr.pos.25.40g = opt.dfr.pos.25[which(opt.dfr.pos.25$p.bm >= 40),]

  plot_ly(opt.dfr.pos.25.40g, x = ~L_nought, y = ~P_nought, z = ~h.bm, type = 'contour')
  
  plot_ly(opt.dfr.pos.25.40g, x = ~L_nought, y = ~P_nought, z = ~profit, type = 'contour')
  plot_ly(opt.dfr.pos.25.40g, x = ~L_nought, y = ~P_nought, z = ~ROI, type = 'contour')

ros.opt = max(opt.dfr.pos.25.40g$ROI)

ros.l0 = opt.dfr.pos.25.40g$L_nought[opt.dfr.pos.25.40g$ROI == ros.opt]  # rosenbergii optimal starting length
ros.p0 = opt.dfr.pos.25.40g$P_nought[opt.dfr.pos.25.40g$ROI == ros.opt]  # rosenbergii optimal starging density
ros.ht = opt.dfr.pos.25.40g$h.t[opt.dfr.pos.25.40g$ROI == ros.opt]       # rosenbergii optimal harvest time
ros.hbm = opt.dfr.pos.25.40g$h.bm[opt.dfr.pos.25.40g$ROI == ros.opt]     # rosenbergii harvest biomass at optimal harvest
ros.pbm = opt.dfr.pos.25.40g$p.bm[opt.dfr.pos.25.40g$ROI == ros.opt]     # rosenbergii prawn mass at optimal harvest


# # Fit density-dependent parameters to data from Ranjeet & Kurup 2010 ##########
# beta.params = function(mu, var) {
#   alpha = ((1 - mu)/var - 1/mu)*mu^2
#   beta = alpha*(1/mu - 1)
#   return(list(alpha = alpha, beta = beta))
# }
# 
# prawn.mle = function(params) {
#   par.aqua['gam'] = params[1]
#   par.aqua['phi'] = params[2]
#   
#   time = seq(0,30*8,1)
#   
#   nstart1 = c(P = 5000, L = 25)
#   output1 = as.data.frame(ode(nstart1, time, prawn_biomass, par.aqua))
#   mean.size1 = ((par.aqua['a']*(output1$L[output1$time == max(time)]/10)^par.aqua['b'])/10)
#   surv1 = output1$P[output1$time == max(time)] / output1$P[output1$time == 0]
#     
#   nstart2 = c(P = 10000, L = 25)
#   output2 = as.data.frame(ode(nstart2, time, prawn_biomass, par.aqua))
#   mean.size2 = ((par.aqua['a']*(output2$L[output2$time == max(time)]/10)^par.aqua['b'])/10)
#   surv2 = output2$P[output2$time == max(time)] / output2$P[output2$time == 0]
#     
#   nstart3 = c(P = 15000, L = 25)
#   output3 = as.data.frame(ode(nstart3, time, prawn_biomass, par.aqua))
#   mean.size3 = ((par.aqua['a']*(output3$L[output3$time == max(time)]/10)^par.aqua['b'])/10)
#   surv3 = output3$P[output3$time == max(time)] / output3$P[output3$time == 0]
#     
#   nstart4 = c(P = 25000, L = 25)
#   output4 = as.data.frame(ode(nstart4, time, prawn_biomass, par.aqua))
#   mean.size4 = ((par.aqua['a']*(output4$L[output4$time == max(time)]/10)^par.aqua['b'])/10)
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


