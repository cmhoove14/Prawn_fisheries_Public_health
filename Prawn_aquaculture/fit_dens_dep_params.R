source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')

t.rk06 = c(1:(8*30)) # time to harvest in Ranjeet and Kurup trials is 8 months 

nsim = 20
surv.df = expand.grid(om = seq(5e-10, 1e-8, length.out = nsim), 
                      gam = seq(5e-7, 1e-5, length.out = nsim),
                      P0 = seq(10000, 70000, length.out = nsim))

surv.df$surv = 0
surv.df$mean_B_perday = 0

par.aqua['k'] = k.ros       # Growth rate (mm/day, M. rosenbergii)

for(k in 1:nrow(surv.df)){
  start = c(P = surv.df[k,3], L = 32)  # set starting density
  par.aqua['gam'] = surv.df[k,2]       # set gamma parameter
  par.aqua['om'] = surv.df[k,1]        # set omega parameter
  
  op = as.data.frame(ode(start,t.rk06,prawn_biomass,par.aqua))   # run simulation
  op$B = par.aqua['a.p'] *op$L^par.aqua['b.p']       # Mean prawn biomass, transformed from length
  
  surv.df[k,4] = op$P[op$time == max(op$time)]/surv.df[k,3]      # store survival rate in simulation
  surv.df[k,5] = op$B[op$time == max(op$time)]/max(op$time)      # store mean daily growth rate in simulation  
  
  if(k %% 1000==0) print(k)
}

for(o in 1:nsim){
  plot(surv.df$P0[surv.df$om == unique(surv.df$om)[o]]/10000, surv.df$surv[surv.df$om == unique(surv.df$om)[o]],
       pch = 16, ylim = c(0,1), xlab = 'Stocking density (P/m^2)', ylab = '%survival',
       main = paste0('omega = ', unique(surv.df$om)[o]))
    points(rk06$dens, rk06$surv, pch = 17, col = 2, cex = 1.25)
    abline(om.lm, lty = 2, col = 2)
}

# 5e-9 fits best
  par.aqua['om'] <- 5e-9

for(g in 1:nsim){
  plot(surv.df$P0[surv.df$gam == unique(surv.df$gam)[g]]/10000, surv.df$mean_B_perday[surv.df$gam == unique(surv.df$gam)[g]],
       pch = 16, ylim = c(0,1), xlab = 'Stocking density (P/m^2)', ylab = 'Mean size',
       main = paste0('gamma = ', unique(surv.df$gam)[g]))
    points(rk06$dens, rk06$mean_B_perday, pch = 17, col = 2, cex = 1.25)
    abline(gam.lm, lty = 2, col = 2)
}

# 3.5e-6 fits best
  par.aqua['gam'] <- 3.5e-6