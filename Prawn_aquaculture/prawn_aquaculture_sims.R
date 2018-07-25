#Generate plot of aquaculture cycle over two years to display in model summary figure
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')

#Simulation for eumetric curve for M. volenhovenii #####
opt.df = expand.grid(L_nought = 33, P_nought = seq(100, 10000, 100)) 
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

    op = as.data.frame(ode(start,t.p,prawn_biomass,par.aqua))
    op$B = ((par.aqua['a.p']*(op$L/10)^par.aqua['b.p'])/10)     # Mean prawn biomass, transformed from length
    
    opt.df[k,3] = op$time[op$B*op$P==max(op$B*op$P)]  # harvest time
    opt.df[k,4] = max(op$B*op$P)                      # Total prawn biomass
    opt.df[k,5] = predict(eta.lm, newdata = data.frame(dens = opt.df[k,2]/area)) # Fraction of harvest that's marketable as function of stocking density
    opt.df[k,6] = opt.df[k,4] * opt.df[k,5]           # Marketable harvest (total biomass times fraction marketable size)
    opt.df[k,7] = op$B[op$B*op$P==max(op$B*op$P)]     # Prawn size at harvest
    opt.df[k,8] = op$P[op$B*op$P==max(op$B*op$P)] / opt.df[k,2]   # Prawn % survival at harvest
    opt.df[k,9] = p*(opt.df[k,4]/1000)*opt.df[k,5]    # Raw Profit
    opt.df[k,10] = p*(opt.df[k,4]/1000)*opt.df[k,5]*exp(-delta*(op$time[op$B*op$P==opt.df[k,4]])) - c*(start["P"])  #Net profit
    opt.df[k,11] = (p*(opt.df[k,4]/1000)*opt.df[k,5]*exp(-delta*(op$time[op$B*op$P==opt.df[k,4]])) - c*(start["P"])) / (c*(start["P"])) # ROI  
    
  }

opt.df$p.t = opt.df$P_nought * opt.df$p.surv     #number of prawns alive at harvest
opt.df$Species = "M. volenhovenii"
  opt.df.criteria = subset(opt.df, p.bm >=30 & h.t <= 365)
  opt.vol = opt.df[which(opt.df$Profit == max(opt.df$Profit)),]

#Simulation for eumetric curve for M. rosenbergii #####
opt.df.ros = expand.grid(L_nought = 33, P_nought = seq(100, 10000, 100)) 
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
    opt.df.ros[k,9] = p*(opt.df.ros[k,4]/1000)*opt.df.ros[k,5]    # Raw Profit
    opt.df.ros[k,10] = p*(opt.df.ros[k,4]/1000)*opt.df.ros[k,5]*exp(-delta*(op$time[op$B*op$P==opt.df.ros[k,4]])) - c*(start["P"])  #Net profit
    opt.df.ros[k,11] = (p*(opt.df.ros[k,4]/1000)*opt.df.ros[k,5]*exp(-delta*(op$time[op$B*op$P==opt.df.ros[k,4]])) - c*(start["P"])) / (c*(start["P"])) # ROI  
    
  }
 
opt.df.ros$p.t = opt.df.ros$P_nought * opt.df.ros$p.surv     #number of prawns alive at harvest
opt.df.ros$Species = "M. rosenbergii"
  opt.df.ros.criteria = subset(opt.df.ros, p.bm >=30 & h.t <= 365)
  opt.ros = opt.df.ros[which(opt.df.ros$Profit == max(opt.df.ros$Profit)),]

    
#Simulate an aquaculture cycle based on profit-optimized stocking density for M. volenhovenii above ###########
t.p = c(0:(365*2))
    
#run model through 2 years for M. volenhovenii ######
par.aqua['k'] = 0.00339726  # alternate value for M. volenhovenii
nstart.p.vol = c(P = opt.vol$P_nought, L = opt.vol$L_nought)   #optimal starting conditions for M. volenhovenii as estimated above
op.vol = as.data.frame(ode(nstart.p.vol,t.p,prawn_biomass,par.aqua))
  c = 0.1+0.005*(nstart.p.vol["L"]-33) # make juvenile prawns cost a function of their size with reference at $0.1/juvenile

#post-process to estimate additional parameters
# 0) prawn density rather than raw prawn numbers
  op.vol$P.dens = op.vol$P / area
# 1) mean prawn biomass (allometric function)
  op.vol$B = ((par.aqua['a.p']*(op.vol$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op.vol$Bt = op.vol$B*op.vol$P / 1000
# 2.5) marketable prawns at harvest
  eta.vol = predict(eta.lm, newdata = data.frame(dens = nstart.p.vol["P"]/area)) # Fraction of harvest that's marketable as function of stocking density
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op.vol$profit = eta.vol*p*(op.vol$Bt/1000)*exp(-delta*(op.vol$t)) - c*(nstart.p["P"])
# 4) starting total biomass
  start.mass.kg.vol = op.vol$Bt[op.vol$time==0]/1000
# 5) harvest mass in kg (harvest assumed to occur when biomass is maximized)   
  harvest.mass.kg.vol = eta.vol*max(op.vol$Bt)/1000
# 6) average mass of prawns at harvest  
  harvest.size.vol = op.vol$B[op.vol$Bt==max(op.vol$Bt)]
# 7) time of harvest (when biomass is maximized)
  harvest.time.vol = op.vol$time[op.vol$Bt==max(op.vol$Bt)]
# 8) time in months
  op.vol$time.mo = op.vol$time / 30
# 9) species
  op.vol$Species = 'M. volenhovenii'

#run model through 2 years for M. rosenbergii ######
par.aqua['k'] = 0.0104333333  # alternate value for M. rosenbergii, from Sampaio & Valenti 1996
nstart.p.ros = c(P = opt.ros$P_nought, L = opt.ros$L_nought)   #optimal starting conditions for M. rosenbergii as estimated above
op.ros = as.data.frame(ode(nstart.p.ros,t.p,prawn_biomass,par.aqua))
  c = 0.1+0.005*(nstart.p.ros["L"]-33) # make juvenile prawns cost a function of their size with reference at $0.1/juvenile

#post-process to estimate additional parameters
# 0) prawn density rather than raw prawn numbers
  op.ros$P.dens = op.ros$P / area
# 1) mean prawn biomass (allometric function)
  op.ros$B = ((par.aqua['a.p']*(op.ros$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op.ros$Bt = op.ros$B*op.ros$P /1000 
# 2.5) marketable prawns at harvest
  eta.ros = predict(eta.lm, newdata = data.frame(dens = nstart.p.ros['P']/area)) # Fraction of harvest that's marketable as function of stocking density
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op.ros$profit = eta.ros*p*(op.ros$Bt/1000)*exp(-delta*(op.ros$t)) - c*(nstart.p["P"]/1000)
# 4) starting total biomass
  start.mass.kg.ros = op.ros$Bt[op.ros$time==0]/1000
# 5) harvest mass in kg (harvest assumed to occur when biomass is maximized)   
  harvest.mass.kg.ros = eta.ros*max(op.ros$Bt)/1000
# 6) average mass of prawns at harvest  
  harvest.size.ros = op.ros$B[op.ros$Bt==max(op.ros$Bt)]
# 7) time of harvest (when biomass is maximized)
  harvest.time.ros = op.ros$time[op.ros$Bt==max(op.ros$Bt)]
# 8) time in months
  op.ros$time.mo = op.ros$time / 30
# 9) species
  op.ros$Species = 'M. rosenbergii'
  
save.image(file = "Prawn_aquaculture/aquaculture_sims.RData")  