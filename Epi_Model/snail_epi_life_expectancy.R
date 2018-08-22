#code used to estimate mean life expectancy in snail population
#original life expectancy among entire snail cohort was 22 days which was causing problems in the integrated 
#model as well as very low snail densities at equalibrium (~4 snails/m^2)
#this code was used to tweak mortality rate parameters to get mean life expectancy to ~40 days

#It was found that snail population dynamics were actually most sensitive to infection dynamics. 
#high transmission led to high mortality rates among the snail population due to accumulation in I compartment
#transmission parameters were tweaked to give mean lfie expectancy of ~44 days

source('Epi_Model/snail_epi_mod_no_diag.R')

nstart.ls = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
              E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
              I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
              Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu)

t2 = c(1:(365*2))

#These parameters were tweaked to increase mean lifespan
#par.snails['muN1'] = 1/50  #Original 1/40
#par.snails['muN2'] = 1/75  #Original 1/50 
#par.snails['muN3'] = 1/100  #Original 1/60
#par.snails['muI'] = 1/25

#rerun the model to equilibrium for accurate estimation of lifespan
sn.eqbm.run = as.data.frame(ode(nstart.ls,t2,snail_epi_allvh,par.snails))
  sn.eqbm = sn.eqbm.run[dim(sn.eqbm.run)[1],]

#plot worm burden trajectory to check it doesn't go crazy  
plot(sn.eqbm.run$time, sn.eqbm.run$Wt, type='l', lwd = 2, col = 2,
     xlab = 'time', ylab = 'mean worm burden')  
  
par.snails['f'] = 0

sn.run.ls = as.data.frame(ode(nstart.ls,t2,snail_epi_allvh,par.snails))
  sn.run.ls$S.t = (sn.run.ls$S1 + sn.run.ls$S2 + sn.run.ls$S3)
  sn.run.ls$E.t = (sn.run.ls$E1 + sn.run.ls$E2 + sn.run.ls$E3)
  sn.run.ls$I.t = (sn.run.ls$I1 + sn.run.ls$I2 + sn.run.ls$I3)
  sn.run.ls$size1 = (sn.run.ls$S1 + sn.run.ls$E1 + sn.run.ls$I1)
  sn.run.ls$size2 = (sn.run.ls$S2 + sn.run.ls$E2 + sn.run.ls$I2)
  sn.run.ls$size3 = (sn.run.ls$S3 + sn.run.ls$E3 + sn.run.ls$I3)
  
  sn.run.ls$N = sn.run.ls$S.t + sn.run.ls$E.t + sn.run.ls$I.t

plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
     xlab = 'time', ylab = 'snail infection dynamics')
  lines(sn.run.ls$time, sn.run.ls$S.t, col = 3, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$E.t, col = 6, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$I.t, col = 2, lwd = 2)
    legend('topright', legend = c('N', 'S', 'E', 'I'), lwd = 2, col = c(1,3,6,2), bty = 'n', cex = 0.7)
 
plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
      xlab = 'time', ylab = 'snail size dynamics')
  lines(sn.run.ls$time, sn.run.ls$size1, col = 7, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$size2, col = 5, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$size3, col = 4, lwd = 2)
    legend('topright', legend = c('N', '1', '2', '3'), lwd = 2, col = c(1,7,5,4), bty = 'n', cex = 0.7)
    
# Calculate PDF f(t) = -dS(t)/dt (i.e. numerically approximate derivative of survival curve)
vect = sn.run.ls$N
deriv = numeric(length(vect)-1)
for (i in 1:(length(vect)-1)) {
  deriv[i] = vect[i] - vect[i+1]
}

deriv2 = deriv/sn.run.ls$N[1]
deriv2 = deriv2[-(length(deriv2))]
sum(deriv2) 

x = 0:(length(deriv2)-1)

plot(x = x, y = deriv2, type = 'l')

# Calculate life expectancy
sum(x*deriv2)    #original ~22 days

#with 1/30, 1/60, 1/90; mean lifespan is ~21 days
#with 1/50, 1/75, 1/100; mean lifespan is ~25 days
#with 1/50, 1/75, 1/100 and muI = 0; mean lifespan is ~71 days

#same procedure, but "turn off" infection
par.snails['beta'] = 0

sn.eqbm.run = as.data.frame(ode(nstart.sn,t.sn,snail_epi_allvh,par.snails))
  sn.eqbm = sn.eqbm.run[dim(sn.eqbm.run)[1],]

#plot worm burden trajectory to check it doesn't go crazy  
plot(sn.eqbm.run$time, sn.eqbm.run$W, type='l', lwd = 2, col = 2,
     xlab = 'time', ylab = 'mean worm burden')  

par.snails['f'] = 0

sn.run.ls = as.data.frame(ode(nstart.ls,t2,snail_epi_allvh,par.snails))
  sn.run.ls$S.t = (sn.run.ls$S1 + sn.run.ls$S2 + sn.run.ls$S3)
  sn.run.ls$E.t = (sn.run.ls$E1 + sn.run.ls$E2 + sn.run.ls$E3)
  sn.run.ls$I.t = (sn.run.ls$I1 + sn.run.ls$I2 + sn.run.ls$I3)
  sn.run.ls$size1 = (sn.run.ls$I1 + sn.run.ls$S1 + sn.run.ls$E1)
  sn.run.ls$size2 = (sn.run.ls$S2 + sn.run.ls$E2 + sn.run.ls$I2)
  sn.run.ls$size3 = (sn.run.ls$S3 + sn.run.ls$E3 + sn.run.ls$I3)

sn.run.ls$N = sn.run.ls$S.t + sn.run.ls$E.t + sn.run.ls$I.t

plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
     xlab = 'time', ylab = 'snail infection dynamics')
  lines(sn.run.ls$time, sn.run.ls$S.t, col = 3, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$E.t, col = 6, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$I.t, col = 2, lwd = 2)
    legend('topright', legend = c('N', 'S', 'E', 'I'), lwd = 2, col = c(1,3,6,2), bty = 'n', cex = 0.7)

plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
     xlab = 'time', ylab = 'snail size dynamics')
  lines(sn.run.ls$time, sn.run.ls$size1, col = 7, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$size2, col = 5, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$size3, col = 4, lwd = 2)
    legend('topright', legend = c('N', '1', '2', '3'), lwd = 2, col = c(1,7,5,4), bty = 'n', cex = 0.7)

# Calculate PDF f(t) = -dS(t)/dt (i.e. numerically approximate derivative of survival curve)
vect = sn.run.ls$N
deriv = numeric(length(vect)-1)
for (i in 1:(length(vect)-1)) {
  deriv[i] = vect[i] - vect[i+1]
}

deriv2 = deriv/sn.run.ls$N[1]
deriv2 = deriv2[-(length(deriv2))]
sum(deriv2) 

x = 0:(length(deriv2)-1)

plot(x = x, y = deriv2, type = 'l')

# Calculate life expectancy
sum(x*deriv2)    #original ~22 days; with infection turned off (beta = 0) lifespan = 42.099
