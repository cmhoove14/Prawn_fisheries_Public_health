## Modeling prawn growth and mortality to maximize biomass

## To-do list: #########
  # Make cost of initial prawn stocking a function of their size (unless we want to fix initial stocking size)
  # Make price of harvested prawns a function of their size such that larger prawns are more valuable
  # Fit density-dependent parameters to data from Ranjeet & Kurup 2010 (or other Macrobrachium fisheries data)
    # Note: may also want to use net production data for fitting?
  # Investigate price-to-cost ratios in other locations for sensitivity analysis on maximum profit estimation

require(plotly)

source('Prawn_aquaculture/prawn_aquaculture_mod.R')
area = 10000 #Working with 1 hectare site
# Run & plot typical volenhovenii cycle from TAMU rosenbergii document ######
  P0.tamu = 2.5*8000      #8000 juveniles stocked per acre, ~2.5 acres per hectare
  L0.tamu = 38            #38mm = ~.78 grams
  c.tamu = 0.07           #7 cents / juvenile; assume that actual juveniles are the only cost
  p.tamu = 5*2.2          #price per kg of adult prawns (TAMU reports $4-$7 per lb)
  
  nstart.tamu = c(P = P0.tamu, L = L0.tamu)
  
  op.tamu = as.data.frame(ode(nstart.tamu,t.p,prawn_biomass,par.aqua))
    op.tamu$B = ((par.aqua['a.p']*(op.tamu$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
    op.tamu$Bt = op.tamu$B*op.tamu$P                                           # Total prawn biomass
    # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
    op.tamu$profit = p.tamu*(op.tamu$Bt/1000)*exp(-delta*(op.tamu$t)) - (c.tamu*nstart.p["P"])   
    #ROI as net profit over cost
    op.tamu$roi = ((p.tamu*(op.tamu$Bt/1000)*exp(-delta*(op.tamu$t))) - 
                    (c.tamu*nstart.p["P"])) /
                  (c.tamu*nstart.p["P"])

      bt0.tamu = op.tamu$Bt[op.tamu$time==0]/1000
      btT.tamu = max(op.tamu$Bt)/1000
      BT.tamu = op.tamu$B[op.tamu$Bt==max(op.tamu$Bt)]  
      PT.tamu = op.tamu$P[op.tamu$Bt==max(op.tamu$Bt)]  
      T.tamu = op.tamu$time[op.tamu$Bt==max(op.tamu$Bt)]
      prof.tamu = op.tamu$profit[op.tamu$Bt==max(op.tamu$Bt)]
      roi.tamu = op.tamu$roi[op.tamu$Bt==max(op.tamu$Bt)]
      
plot(x = op.tamu$time, y = op.tamu$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(op.tamu$time)), ylim = c(0, max(op.tamu$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.tamu[2]), ' mm)', sep = ''))
  lines(x = op.tamu$time, y = op.tamu$Bt/1000, col = 'blue', lwd=2)
  lines(x = op.tamu$time, y = op.tamu$B, col = 'purple', lwd=2, lty=2)
  lines(x = op.tamu$time, y = op.tamu$L, col = 'green', lwd=2)
  abline(v = T.tamu, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5, bty = 'n')
  legend('bottomright', legend = c(paste('Starting mass =', round(bt0.tamu), 'kg', sep = ' '),
                           paste('Total harvest mass =', round(btT.tamu), 'kg', sep = ' '), 
                           paste('Mean harvest mass =', round(BT.tamu), 'g', sep = ' '), 
                           paste('Time of harvest =', round(T.tamu), 'days', sep = ' ')), cex=0.5, bty = 'n')
  
plot(op.tamu$time, op.tamu$profit, type = 'l', col = 'darkgreen', lwd = 2,
     xlab = 'time', ylab = 'profit')  
plot(op.tamu$time, op.tamu$roi, type = 'l', col = 'blue', lwd = 2,
     xlab = 'time', ylab = 'ROI')

# Run & plot typical rosenbergii cycle from TAMU Rosenbergii document ######
  P0.tamu = 2.5*8000      #8000 juveniles stocked per acre, ~2.5 acres per hectare
  L0.tamu = 38            #38mm = ~.78 grams
  c.tamu = 0.07           #7 cents / juvenile; assume that actual juveniles are the only cost
  p.tamu = 5*2.2          #price per kg of adult prawns (TAMU reports $4-$7 per lb)
  
  par.aqua['k'] = 0.0104333333   #alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 
  nstart.tamu = c(P = P0.tamu, L = L0.tamu)
  
  op.tamu.ros = as.data.frame(ode(nstart.tamu,t.p,prawn_biomass,par.aqua))
  op.tamu.ros$B = ((par.aqua['a.p']*(op.tamu.ros$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
  op.tamu.ros$Bt = op.tamu.ros$B*op.tamu.ros$P                                           # Total prawn biomass
  # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
  op.tamu.ros$profit = p.tamu*(op.tamu.ros$Bt/1000)*exp(-delta*(op.tamu.ros$t)) - (c.tamu*nstart.p["P"])   
  #ROI as net profit over cost
  op.tamu.ros$roi = ((p.tamu*(op.tamu.ros$Bt/1000)*exp(-delta*(op.tamu.ros$t))) - 
                   (c.tamu*nstart.p["P"])) /
                    (c.tamu*nstart.p["P"])

    bt0.tamu.ros = op.tamu.ros$Bt[op.tamu.ros$time==0]/1000
    btT.tamu.ros = max(op.tamu.ros$Bt)/1000
    BT.tamu.ros = op.tamu.ros$B[op.tamu.ros$Bt==max(op.tamu.ros$Bt)]  
    PT.tamu.ros = op.tamu.ros$P[op.tamu.ros$Bt==max(op.tamu.ros$Bt)]  
    T.tamu.ros = op.tamu.ros$time[op.tamu.ros$Bt==max(op.tamu.ros$Bt)]
    prof.tamu.ros = op.tamu.ros$profit[op.tamu.ros$Bt==max(op.tamu.ros$Bt)]
    roi.tamu.ros = op.tamu.ros$roi[op.tamu.ros$Bt==max(op.tamu.ros$Bt)]

plot(x = op.tamu.ros$time, y = op.tamu.ros$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(op.tamu.ros$time)), ylim = c(0, max(op.tamu.ros$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.tamu.ros[2]), ' mm)', sep = ''))
  lines(x = op.tamu.ros$time, y = op.tamu.ros$Bt/1000, col = 'blue', lwd=2)
  lines(x = op.tamu.ros$time, y = op.tamu.ros$B, col = 'purple', lwd=2, lty=2)
  lines(x = op.tamu.ros$time, y = op.tamu.ros$L, col = 'green', lwd=2)
  abline(v = T.tamu.ros, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5, bty = 'n')
  legend('bottomright', legend = c(paste('Starting mass =', round(bt0.tamu.ros), 'kg', sep = ' '),
                                   paste('Total harvest mass =', round(btT.tamu.ros), 'kg', sep = ' '), 
                                   paste('Mean harvest mass =', round(BT.tamu.ros), 'g', sep = ' '), 
                                   paste('Time of harvest =', round(T.tamu.ros), 'days', sep = ' ')), cex=0.5, bty = 'n')

plot(op.tamu.ros$time, op.tamu.ros$profit, type = 'l', col = 'darkgreen', lwd = 2,
     xlab = 'time', ylab = 'profit')  
plot(op.tamu.ros$time, op.tamu.ros$roi, type = 'l', col = 'blue', lwd = 2,
     xlab = 'time', ylab = 'ROI')

# Run & plot cycles with different stocking parameters from Willis and Berrigan 1977 paper, scenario 1############
  P0.wb1 = 5*area      #5 juveniles per square meter
  L0.wb1 = 38          #38mm = ~.78 grams
  c.wb1 = 0.07         #7 cents / juvenile; assume that actual juveniles are the only cost
  p.wb1= 5*2.2        #price per kg of adult prawns (TAMU reports $4-$7 per lb)
  
  nstart.wb1 = c(P = P0.wb1, L = L0.wb1)

op.wb1 = as.data.frame(ode(nstart.wb1,t.p,prawn_biomass,par.aqua))
  op.wb1$B = ((par.aqua['a.p']*(op.wb1$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
  op.wb1$Bt = op.wb1$B*op.wb1$P                                           # Total prawn biomass
# Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
  op.wb1$profit = p.wb1*(op.wb1$Bt/1000)*exp(-delta*(op.wb1$t)) - (c.wb1*nstart.p["P"])   
# ROI as net profit over cost
  op.wb1$roi = ((p.wb1*(op.wb1$Bt/1000)*exp(-delta*(op.wb1$t))) - 
                   (c.wb1*nstart.p["P"])) /
                (c.wb1*nstart.p["P"])
  
  bt0.wb1 = op.wb1$Bt[op.wb1$time==0]/1000
  T.wb1 = 168
  btT.wb1 = op.wb1$Bt[op.wb1$time == T.wb1]/1000
  BT.wb1 = op.wb1$B[op.wb1$time == T.wb1]  
  PT.wb1 = op.wb1$P[op.wb1$time == T.wb1]  
  prof.wb1 = op.wb1$profit[op.wb1$time == T.wb1]
  roi.wb1 = op.wb1$roi[op.wb1$time == T.wb1]
  
plot(x = op.wb1$time, y = op.wb1$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(op.wb1$time)), ylim = c(0, max(op.wb1$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.wb1[2]), ' mm)', sep = ''))
  lines(x = op.wb1$time, y = op.wb1$Bt/1000, col = 'blue', lwd=2)
  lines(x = op.wb1$time, y = op.wb1$B, col = 'purple', lwd=2, lty=2)
  lines(x = op.wb1$time, y = op.wb1$L, col = 'green', lwd=2)
  abline(v = T.wb1, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5, bty = 'n')
  legend('bottomright', legend = c(paste('Starting mass =', round(bt0.wb1), 'kg', sep = ' '),
                                   paste('Total harvest mass =', round(btT.wb1), 'kg', sep = ' '), 
                                   paste('Mean harvest mass =', round(BT.wb1), 'g', sep = ' '), 
                                   paste('Time of harvest =', round(T.wb1), 'days', sep = ' ')), cex=0.5, bty = 'n')
  
  plot(op.wb1$time, op.wb1$profit, type = 'l', col = 'darkgreen', lwd = 2,
       xlab = 'time', ylab = 'profit')  
  plot(op.wb1$time, op.wb1$roi, type = 'l', col = 'blue', lwd = 2,
       xlab = 'time', ylab = 'ROI')
  
# Run & plot cycles with different stocking parameters from Willis and Berrigan 1977 paper, scenario 2############
  P0.wb2 = 5*area      #5 post larvae per square meter
  L0.wb2 = 17          #17mm = ~.056 grams
  c.wb2 = 0.007       #0.7 cents / postlarvae; assume that actual postlarvae are the only cost; assume post-larvae are 10% cost of juveniles
  p.wb2= 5*2.2        #price per kg of adult prawns (TAMU reports $4-$7 per lb)
  
  nstart.wb2 = c(P = P0.wb2, L = L0.wb2)
  
  op.wb2 = as.data.frame(ode(nstart.wb2,t.p,prawn_biomass,par.aqua))
  op.wb2$B = ((par.aqua['a.p']*(op.wb2$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
  op.wb2$Bt = op.wb2$B*op.wb2$P                                           # Total prawn biomass
  # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
  op.wb2$profit = p.wb2*(op.wb2$Bt/1000)*exp(-delta*(op.wb2$t)) - (c.wb2*nstart.p["P"])   
  # ROI as net profit over cost
  op.wb2$roi = ((p.wb2*(op.wb2$Bt/1000)*exp(-delta*(op.wb2$t))) - 
                  (c.wb2*nstart.p["P"])) /
    (c.wb2*nstart.p["P"])
  
  bt0.wb2 = op.wb2$Bt[op.wb2$time==0]/1000
  T.wb2 = 168
  btT.wb2 = op.wb2$Bt[op.wb2$time == T.wb2]/1000
  BT.wb2 = op.wb2$B[op.wb2$time == T.wb2]  
  PT.wb2 = op.wb2$P[op.wb2$time == T.wb2]  
  prof.wb2 = op.wb2$profit[op.wb2$time == T.wb2]
  roi.wb2 = op.wb2$roi[op.wb2$time == T.wb2]
  
plot(x = op.wb2$time, y = op.wb2$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(op.wb2$time)), ylim = c(0, max(op.wb2$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.wb2[2]), ' mm)', sep = ''))
  lines(x = op.wb2$time, y = op.wb2$Bt/1000, col = 'blue', lwd=2)
  lines(x = op.wb2$time, y = op.wb2$B, col = 'purple', lwd=2, lty=2)
  lines(x = op.wb2$time, y = op.wb2$L, col = 'green', lwd=2)
  abline(v = T.wb2, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5, bty = 'n')
  legend('bottomright', legend = c(paste('Starting mass =', round(bt0.wb2), 'kg', sep = ' '),
                                   paste('Total harvest mass =', round(btT.wb2), 'kg', sep = ' '), 
                                   paste('Mean harvest mass =', round(BT.wb2), 'g', sep = ' '), 
                                   paste('Time of harvest =', round(T.wb2), 'days', sep = ' ')), cex=0.5, bty = 'n')
  
  plot(op.wb2$time, op.wb2$profit, type = 'l', col = 'darkgreen', lwd = 2,
       xlab = 'time', ylab = 'profit')  
  plot(op.wb2$time, op.wb2$roi, type = 'l', col = 'blue', lwd = 2,
       xlab = 'time', ylab = 'ROI')
  
  
# Run & plot cycles with different stocking parameters from Willis and Berrigan 1977 paper, scenario 3############
  P0.wb3 = 10*area      #5 post larvae per square meter
  L0.wb3 = 17          #17mm = ~.056 grams
  c.wb3 = 0.007       #0.7 cents / postlarvae; assume that actual postlarvae are the only cost; assume post-larvae are 10% cost of juveniles
  p.wb3= 5*2.2        #price per kg of adult prawns (TAMU reports $4-$7 per lb)
  
  nstart.wb3 = c(P = P0.wb3, L = L0.wb3)
  
op.wb3 = as.data.frame(ode(nstart.wb3,t.p,prawn_biomass,par.aqua))
  op.wb3$B = ((par.aqua['a.p']*(op.wb3$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
  op.wb3$Bt = op.wb3$B*op.wb3$P                                           # Total prawn biomass
  # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
  op.wb3$profit = p.wb3*(op.wb3$Bt/1000)*exp(-delta*(op.wb3$t)) - (c.wb3*nstart.p["P"])   
  # ROI as net profit over cost
  op.wb3$roi = ((p.wb3*(op.wb3$Bt/1000)*exp(-delta*(op.wb3$t))) - 
                  (c.wb3*nstart.p["P"])) /
    (c.wb3*nstart.p["P"])
  
  bt0.wb3 = op.wb3$Bt[op.wb3$time==0]/1000
  T.wb3 = 168
  btT.wb3 = op.wb3$Bt[op.wb3$time == T.wb3]/1000
  BT.wb3 = op.wb3$B[op.wb3$time == T.wb3]  
  PT.wb3 = op.wb3$P[op.wb3$time == T.wb3]  
  prof.wb3 = op.wb3$profit[op.wb3$time == T.wb3]
  roi.wb3 = op.wb3$roi[op.wb3$time == T.wb3]
  
plot(x = op.wb3$time, y = op.wb3$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(op.wb3$time)), ylim = c(0, max(op.wb3$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.wb3[2]), ' mm)', sep = ''))
  lines(x = op.wb3$time, y = op.wb3$Bt/1000, col = 'blue', lwd=2)
  lines(x = op.wb3$time, y = op.wb3$B, col = 'purple', lwd=2, lty=2)
  lines(x = op.wb3$time, y = op.wb3$L, col = 'green', lwd=2)
  abline(v = T.wb3, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5, bty = 'n')
  legend('bottomright', legend = c(paste('Starting mass =', round(bt0.wb3), 'kg', sep = ' '),
                                   paste('Total harvest mass =', round(btT.wb3), 'kg', sep = ' '), 
                                   paste('Mean harvest mass =', round(BT.wb3), 'g', sep = ' '), 
                                   paste('Time of harvest =', round(T.wb3), 'days', sep = ' ')), cex=0.5, bty = 'n')
  
  plot(op.wb3$time, op.wb3$profit, type = 'l', col = 'darkgreen', lwd = 2,
       xlab = 'time', ylab = 'profit')  
  plot(op.wb3$time, op.wb3$roi, type = 'l', col = 'blue', lwd = 2,
       xlab = 'time', ylab = 'ROI')
  
  
# Run & plot cycles with different stocking parameters from Willis and Berrigan 1977 paper, scenario 4############
  P0.wb4 = 20*area      #5 post larvae per square meter
  L0.wb4 = 17          #17mm = ~.056 grams
  c.wb4 = 0.007       #0.7 cents / postlarvae; assume that actual postlarvae are the only cost; assume post-larvae are 10% cost of juveniles
  p.wb4= 5*2.2        #price per kg of adult prawns (TAMU reports $4-$7 per lb)
  
  nstart.wb4 = c(P = P0.wb4, L = L0.wb4)
  
op.wb4 = as.data.frame(ode(nstart.wb4,t.p,prawn_biomass,par.aqua))
  op.wb4$B = ((par.aqua['a.p']*(op.wb4$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
  op.wb4$Bt = op.wb4$B*op.wb4$P                                           # Total prawn biomass
  # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
  op.wb4$profit = p.wb4*(op.wb4$Bt/1000)*exp(-delta*(op.wb4$t)) - (c.wb4*nstart.p["P"])   
  # ROI as net profit over cost
  op.wb4$roi = ((p.wb4*(op.wb4$Bt/1000)*exp(-delta*(op.wb4$t))) - 
                  (c.wb4*nstart.p["P"])) /
    (c.wb4*nstart.p["P"])
  
  bt0.wb4 = op.wb4$Bt[op.wb4$time==0]/1000
  T.wb4 = 168
  btT.wb4 = op.wb4$Bt[op.wb4$time == T.wb4]/1000
  BT.wb4 = op.wb4$B[op.wb4$time == T.wb4]  
  PT.wb4 = op.wb4$P[op.wb4$time == T.wb4]  
  prof.wb4 = op.wb4$profit[op.wb4$time == T.wb4]
  roi.wb4 = op.wb4$roi[op.wb4$time == T.wb4]
  
plot(x = op.wb4$time, y = op.wb4$P/100, col = 'red', xlab = 'Time (days)', ylab = 'State variables', 
     type = 'l', lwd=2, xlim = c(0, max(op.wb4$time)), ylim = c(0, max(op.wb4$Bt/1000)+50),
     main = paste('Prawn fishery dynamics\n', '(mean start size = ', as.numeric(nstart.wb4[2]), ' mm)', sep = ''))
  lines(x = op.wb4$time, y = op.wb4$Bt/1000, col = 'blue', lwd=2)
  lines(x = op.wb4$time, y = op.wb4$B, col = 'purple', lwd=2, lty=2)
  lines(x = op.wb4$time, y = op.wb4$L, col = 'green', lwd=2)
  abline(v = T.wb4, lty = 2, lwd = 2)
  legend('topright', legend = c('Prawns (100s)', 'Total biomass (kg)', 'Mean size (g)', 'Mean length (mm)'), 
         lty = 1, col = c('red', 'blue', 'purple', 'green'), cex = 0.5, bty = 'n')
  legend('bottomright', legend = c(paste('Starting mass =', round(bt0.wb4), 'kg', sep = ' '),
                                   paste('Total harvest mass =', round(btT.wb4), 'kg', sep = ' '), 
                                   paste('Mean harvest mass =', round(BT.wb4), 'g', sep = ' '), 
                                   paste('Time of harvest =', round(T.wb4), 'days', sep = ' ')), cex=0.5, bty = 'n')
  
  plot(op.wb4$time, op.wb4$profit, type = 'l', col = 'darkgreen', lwd = 2,
       xlab = 'time', ylab = 'profit')  
  plot(op.wb4$time, op.wb4$roi, type = 'l', col = 'blue', lwd = 2,
       xlab = 'time', ylab = 'ROI')
  
  
# Optimize stocking density based on stocking juveniles ########
  op.juv.dens = seq(.1, 50, length.out = 100)
  
  T.juv = numeric(length(op.juv.dens))
  bt0.juv = numeric(length(op.juv.dens))
  btT.juv = numeric(length(op.juv.dens))
  BT.juv = numeric(length(op.juv.dens)) 
  PT.juv = numeric(length(op.juv.dens)) 
  prof.juv = numeric(length(op.juv.dens))
  roi.juv = numeric(length(op.juv.dens))
  
  c.juv = c.tamu
  p.juv = p.tamu
  
  for(j in 1:length(op.juv.dens)){
    nstart.op.juv = c(P = op.juv.dens[j]*area, L = 38)
    
  op.juv = as.data.frame(ode(nstart.op.juv,t.p,prawn_biomass,par.aqua))
    op.juv$B = ((par.aqua['a.p']*(op.juv$L/10)^par.aqua['b.p'])/10)         # Mean prawn biomass, transformed from length
    op.juv$Bt = op.juv$B*op.juv$P                                           # Total prawn biomass
  # Profit function, in terms of revenue (discounted by time since stocking) minus stocking costs 
    op.juv$profit = p.juv*(op.juv$Bt/1000)*exp(-delta*(op.juv$t)) - (c.juv*nstart.p["P"])   
  # ROI as net profit over cost
    op.juv$roi = ((p.juv*(op.juv$Bt/1000)*exp(-delta*(op.juv$t))) - 
                    (c.juv*nstart.p["P"])) /
                  (c.juv*nstart.p["P"])
    T.this = op.juv$time[op.juv$Bt == max(op.juv$Bt)]
    T.juv[j] = T.this
    bt0.juv[j] = op.juv$Bt[op.juv$time==0]/1000
    btT.juv[j] = op.juv$Bt[op.juv$time == T.this]/1000
    BT.juv[j] = op.juv$B[op.juv$time == T.this]  
    PT.juv[j] = op.juv$P[op.juv$time == T.this]  
    prof.juv[j] = op.juv$profit[op.juv$time == T.this]
    roi.juv[j] = op.juv$roi[op.juv$time == T.this]
  }
  
  plot(op.juv.dens, T.juv, type = 'l', lwd = 2,
       xlab = 'Stocking density', ylab = 'Harvest time')
  
  plot(op.juv.dens, BT.juv, type = 'l', lwd = 2,
       xlab = 'Stocking density', ylab = 'Harvest prawn size')
  
  plot(op.juv.dens, btT.juv, type = 'l', lwd = 2,
       xlab = 'Stocking density', ylab = 'Harvest yield')
  
  plot(op.juv.dens, prof.juv, type = 'l', lwd = 2,
       xlab = 'Stocking density', ylab = 'Profit')
  
  plot(op.juv.dens, roi.juv, type = 'l', lwd = 2,
       xlab = 'Stocking density', ylab = 'ROI')
  
# Plot amount of harvest, biomass multiplier, optimal harvest time against initial stocking density #######
x = c(1:50)
harvest = numeric(length(x))
harvest.m = numeric(length(x))
harvest.t = numeric(length(x))
for (i in x) {
  nstart.d = c(P = 1000*i, L = 37)
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
    