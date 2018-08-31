#code used to estimate mean life expectancy in snail population
#original life expectancy among entire snail cohort was 22 days which was causing problems in the integrated 
#model as well as very low snail densities at equalibrium (~4 snails/m^2)
#this code was used to tweak mortality rate parameters to get mean life expectancy to ~40 days

#It was found that snail population dynamics were actually most sensitive to infection dynamics. 
#high transmission led to high mortality rates among the snail population due to accumulation in I compartment
#transmission parameters were tweaked to give mean lfie expectancy of ~44 days

source('Data/snail_epi_parameters.R')
source("R/Epi_Model/snail_epi_mod_no_diag_immigration.R")

require(tidyverse)

eqbm.run = as.data.frame(ode(nstart.sn,
                             t.sn,
                             snail_epi_allvh_imm,
                             par.snails.imm))

  eqbm = eqbm.run %>% filter(time == max(time)) %>% select(S1:Wu) %>% unlist()

t2 = c(1:(365*2))

#These parameters were tweaked to increase mean lifespan
#par.snails['muN1'] = 1/50  #Original 1/40
#par.snails['muN2'] = 1/75  #Original 1/50 
#par.snails['muN3'] = 1/100  #Original 1/60
#par.snails['muI'] = 1/25

#rerun the model to equilibrium for accurate estimation of lifespan
sn.eqbm.run = as.data.frame(ode(eqbm,
                                t2,
                                snail_epi_allvh_imm,
                                par.snails.imm))

  sn.eqbm = sn.eqbm.run[dim(sn.eqbm.run)[1],]

#plot worm burden trajectory to check it doesn't go crazy  
plot(sn.eqbm.run$time, sn.eqbm.run$Wt, type='l', lwd = 2, col = 6,
     xlab = 'time', ylab = 'mean worm burden')  
  
par.snails.imm['f'] = 0

sn.run.ls = as.data.frame(ode(eqbm,
                              t2,
                              snail_epi_allvh_imm,
                              par.snails.imm)) %>% 
    mutate(S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3  

plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
     xlab = 'time', ylab = 'snail infection dynamics')
  lines(sn.run.ls$time, sn.run.ls$S.t, col = 3, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$E.t, col = 6, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$I.t, col = 2, lwd = 2)
    legend('topright', legend = c('N', 'S', 'E', 'I'), lwd = 2, col = c(1,3,6,2), bty = 'n', cex = 0.7)
 
plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
      xlab = 'time', ylab = 'snail size dynamics')
  lines(sn.run.ls$time, sn.run.ls$t.1, col = 7, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$t.2, col = 5, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$t.3, col = 4, lwd = 2)
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
#with 1/50, 1/75, 1/100; mean lifespan is ~40 days

#same procedure, but "turn off" infection ###########
par.snails.imm['beta'] = 0

#rerun the model to equilibrium for accurate estimation of lifespan
sn.eqbm.run = as.data.frame(ode(eqbm,
                                t2,
                                snail_epi_allvh_imm,
                                par.snails.imm))

  sn.eqbm = sn.eqbm.run[dim(sn.eqbm.run)[1],]

#plot worm burden trajectory to check it doesn't go crazy  
plot(sn.eqbm.run$time, sn.eqbm.run$Wt, type='l', lwd = 2, col = 6,
     xlab = 'time', ylab = 'mean worm burden')  
  
par.snails.imm['f'] = 0

sn.run.ls = as.data.frame(ode(eqbm,
                              t2,
                              snail_epi_allvh_imm,
                              par.snails.imm)) %>% 
    mutate(S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3  

plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
     xlab = 'time', ylab = 'snail infection dynamics')
  lines(sn.run.ls$time, sn.run.ls$S.t, col = 3, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$E.t, col = 6, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$I.t, col = 2, lwd = 2)
    legend('topright', legend = c('N', 'S', 'E', 'I'), lwd = 2, col = c(1,3,6,2), bty = 'n', cex = 0.7)
 
plot(sn.run.ls$time, sn.run.ls$N, type = 'l', lwd = 2,
      xlab = 'time', ylab = 'snail size dynamics')
  lines(sn.run.ls$time, sn.run.ls$t.1, col = 7, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$t.2, col = 5, lwd = 2)
  lines(sn.run.ls$time, sn.run.ls$t.3, col = 4, lwd = 2)
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
