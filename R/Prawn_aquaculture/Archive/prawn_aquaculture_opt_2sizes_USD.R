## Modeling prawn growth and mortality to maximize biomass

## To-do list: #########
  # Make cost of initial prawn stocking a function of their size (unless we want to fix initial stocking size)
  # Make price of harvested prawns a function of their size such that larger prawns are more valuable
  # Fit density-dependent parameters to data from Ranjeet & Kurup 2010 (or other Macrobrachium fisheries data)
    # Note: may also want to use net production data for fitting?
  # Investigate price-to-cost ratios in other locations for sensitivity analysis on maximum profit estimation

require(plotly)

source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')

area = 1000 #Working with 1 sq km site
# Simulate across stocking densities to estimate optimal P_nought for M. volenhovenii ############
opt.df = expand.grid(L_nought = c(33, 67), P_nought = seq(1000, 10000, 500)) 
  opt.df$h.t = 0
  opt.df$h.bm = 0
  opt.df$h.frac = 0
  opt.df$m.bm = 0
  opt.df$p.bm = 0
  opt.df$p.surv = 0
  opt.df$Revenue = 0
  opt.df$Profit = 0
  opt.df$ROI = 0

par.aqua['k'] = 0.00339726       # Growth rate (mm/day), from Nwosu & Wolfi 2006 (M. vollenhovenii)

  for(k in 1:nrow(opt.df)){
    start = c(P = opt.df[k,2], L = opt.df[k,1])
    c = 0.1+0.005*(start["L"]-38) # make juvenile prawns cost a function of their size with reference at $0.1/juvenile
    
    op = as.data.frame(ode(start,t.p,prawn_biomass,par.aqua))
    op$B = ((par.aqua['a.p']*(op$L/10)^par.aqua['b.p'])/10)     # Mean prawn biomass, transformed from length
    
    opt.df[k,3] = op$time[op$B*op$P==max(op$B*op$P)]  # harvest time
    opt.df[k,4] = max(op$B*op$P)                      # Total prawn biomass
    opt.df[k,5] = predict(eta.lm, newdata = data.frame(dens = opt.df[k,2]/area)) # Fraction of harvest that's marketable as function of stocking density
    opt.df[k,6] = opt.df[k,4] * opt.df[k,5]           # Marketable harvest (total biomass times fraction marketable size)
    opt.df[k,7] = op$B[op$B*op$P==max(op$B*op$P)]     # Prawn size at harvest
    opt.df[k,8] = op$P[op$B*op$P==max(op$B*op$P)] / opt.df[k,2]   # Prawn % survival at harvest
    opt.df[k,9] = p*(max(op$B*op$P)/1000)    # Raw Profit
    opt.df[k,10] = p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"])  #Net profit
    opt.df[k,11] = (p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"])) / (c*(start["P"])) # ROI  
    
    if(k %% 100==0) print(k)
  }

opt.df$p.t = opt.df$P_nought * opt.df$p.surv     #number of prawns alive at harvest
opt.df$Size = rep(c('Juvenile', 'Small Adult'), nrow(opt.df)/2)
  opt.df.criteria = subset(opt.df, p.bm >=30 & h.t <= 365)
  opt.vol = opt.df[which(opt.df$Profit == max(opt.df$Profit)),]

#harvest time for two different stocking sizes #########
ggplot(opt.df, aes(x = P_nought)) + 
  theme_bw() +
  labs(x = 'stocking density (P/m^2)', y = 'Harvest time',
       title = 'M. vollenhovenii') +
  geom_hline(yintercept = 365) +
  geom_line(aes(y = h.t, col = Size), size = 1.25) 

#size at harvest for two different stocking sizes #########
ggplot(opt.df, aes(x = P_nought)) + 
  theme_bw() +
  labs(x = 'stocking density (P/m^2)', y = 'Harvest prawn size (g)',
       title = 'M. vollenhovenii') +
  geom_hline(yintercept = 30) +
  geom_line(aes(y = p.bm, col = Size), size = 1.25) 

#harvest biomass for two different stocking sizes #########
ggplot(opt.df, aes(x = P_nought)) + 
  theme_bw() +
  labs(x = 'stocking density (P/m^2)', y = 'Harvest biomass (kg)',
       title = 'M. vollenhovenii') +
  #geom_hline(yintercept = 30) +
  geom_line(aes(y = m.bm/1000, col = Size), size = 1.25) 

#ROI for two different stocking sizes ##########
ggplot(opt.df, aes(x = P_nought)) + 
  theme_bw() +
  labs(x = 'stocking density (P/m^2)', y = 'ROI',
       title = 'M. vollenhovenii') +
  geom_hline(yintercept = 1) +
  geom_line(aes(y = ROI, col = Size), size = 1.25) 

#profit for two different stocking sizes ##########
  ggplot(opt.df, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'Profit (USD)',
         title = 'M. vollenhovenii') +
    geom_hline(yintercept = 1) +
    geom_line(aes(y = Profit, col = Size), size = 1.25) 
  
  
  vol.l0 = opt.vol$L_nought  # volenhovenii optimal starting length
  vol.p0 = opt.vol$P_nought  # volenhovenii optimal starging density
  vol.ht = opt.vol$h.t       # volenhovenii optimal harvest time
  vol.mbm = opt.vol$m.bm/1000# volenhovenii marketable biomass at optimal harvest
  vol.pbm = opt.vol$p.bm     # volenhovenii prawn mass at optimal harvest
  vol.roi = opt.vol$ROI      # volenhovenii ROI at harvest
  
###########  
###########  
# Simulate across stocking densities to estimate optimal P_nought for M. rosenbergii ############
opt.df.ros = expand.grid(L_nought = c(33, 67), P_nought = seq(1000, 10000, 500)) 
  opt.df.ros$h.t = 0
  opt.df.ros$h.bm = 0
  opt.df.ros$h.frac = 0
  opt.df.ros$m.bm = 0
  opt.df.ros$p.bm = 0
  opt.df.ros$p.surv = 0
  opt.df.ros$Revenue = 0
  opt.df.ros$Profit = 0
  opt.df.ros$ROI = 0
  
  par.aqua['k'] = 0.0104333333    # alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333
  
  for(k in 1:nrow(opt.df.ros)){
    start = c(P = opt.df.ros[k,2], L = opt.df.ros[k,1])
    c = 0.1+0.005*(start["L"]-38) # make juvenile prawns cost a function of their size with reference at $0.1/juvenile
    
    op = as.data.frame(ode(start,t.p,prawn_biomass,par.aqua))
    op$B = ((par.aqua['a.p']*(op$L/10)^par.aqua['b.p'])/10)     # Mean prawn biomass, transformed from length
    
    opt.df.ros[k,3] = op$time[op$B*op$P==max(op$B*op$P)]  # harvest time
    opt.df.ros[k,4] = max(op$B*op$P)                      # Total prawn biomass
    opt.df.ros[k,5] = predict(eta.lm, newdata = data.frame(dens = opt.df.ros[k,2]/area)) # Fraction of harvest that's marketable as function of stocking density
    opt.df.ros[k,6] = opt.df.ros[k,4] * opt.df.ros[k,5]           # Marketable harvest (total biomass times fraction marketable size)
    opt.df.ros[k,7] = op$B[op$B*op$P==max(op$B*op$P)]     # Prawn size at harvest
    opt.df.ros[k,8] = op$P[op$B*op$P==max(op$B*op$P)] / opt.df.ros[k,2]   # Prawn % survival at harvest
    opt.df.ros[k,9] = p*(max(op$B*op$P)/1000)    # Raw Profit
    opt.df.ros[k,10] = p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"])  #Net profit
    opt.df.ros[k,11] = (p*(max(op$B*op$P)/1000)*exp(-delta*(op$time[op$B*op$P==max(op$B*op$P)])) - c*(start["P"])) / (c*(start["P"])) # ROI  
    
    if(k %% 100==0) print(k)
  }
 
opt.df.ros$p.t = opt.df.ros$P_nought * opt.df.ros$p.surv     #number of prawns alive at harvest
opt.df.ros$Size = rep(c('Juvenile', 'Small Adult'), nrow(opt.df)/2)
  opt.df.ros.criteria = subset(opt.df.ros, p.bm >=30 & h.t <= 365)
  opt.ros = opt.df.ros[which(opt.df.ros$Profit == max(opt.df.ros$Profit)),]

#harvest time for two different stocking sizes #########
  ggplot(opt.df.ros, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'Harvest time',
         title = 'M. rosenbergii') +
    geom_hline(yintercept = 365) +
    geom_line(aes(y = h.t, col = Size), size = 1.25)
  
#size at harvest for two different stocking sizes #########
  ggplot(opt.df.ros, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'Harvest prawn size (g)',
         title = 'M. rosenbergii') +
    geom_hline(yintercept = 30) +
    geom_line(aes(y = p.bm, col = Size), size = 1.25) 

#harvest biomass for two different stocking sizes #########
  ggplot(opt.df.ros, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'Harvest biomass (kg)',
         title = 'M. rosenbergii') +
    #geom_hline(yintercept = 30) +
    geom_line(aes(y = m.bm/1000, col = Size), size = 1.25) 
 
#ROI for two different stocking sizes ##########
  ggplot(opt.df.ros, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'ROI',
         title = 'M. rosenbergii') +
    geom_hline(yintercept = 1) +
    geom_line(aes(y = ROI, col = Size), size = 1.25) 
 
#Profit for two different stocking sizes ##########
  ggplot(opt.df.ros, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'Profit (USD)',
         title = 'M. rosenbergii') +
    geom_hline(yintercept = 1) +
    geom_line(aes(y = Profit, col = Size), size = 1.25) 

#Revenue for two different stocking sizes ##########
  ggplot(opt.df.ros, aes(x = P_nought)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'Revenue (USD)',
         title = 'M. rosenbergii') +
    geom_hline(yintercept = 1) +
    geom_line(aes(y = Revenue, col = Size), size = 1.25) 
  
