## Epidemiological model including snail size classes; no predation

## To-do list:
  # Continue tuning infection parameters, draw estimates from literature as necessary
  # May want to consider a few aspects of snail dynamics further, e.g.: 
      # Should exposed snails be able to reproduce at a reduced rate?
      # Should growth rates differ by infection status?
      # Think a little more about diagonal transitions
      # Reconsider area scaling
  # Come up with R0 expression?

source('Epi_Model/snail_epi_mod_no_diag.R')

# Run & plot #########
  allvh.eqbm.run$S.t = (allvh.eqbm.run$S1 + allvh.eqbm.run$S2 + allvh.eqbm.run$S3) / area        # density susceptible snails
  allvh.eqbm.run$E.t = (allvh.eqbm.run$E1 + allvh.eqbm.run$E2 + allvh.eqbm.run$E3) / area        # density exposed snails 
  allvh.eqbm.run$I.t = (allvh.eqbm.run$I1 + allvh.eqbm.run$I2 + allvh.eqbm.run$I3 ) / area                         # density infected snails
  allvh.eqbm.run$N.t = (allvh.eqbm.run$S.t + allvh.eqbm.run$E.t + allvh.eqbm.run$I.t)            # density snails
  allvh.eqbm.run$t.1 = (allvh.eqbm.run$S1 + allvh.eqbm.run$E1) / area                          # density snails of size class 1
  allvh.eqbm.run$t.2 = (allvh.eqbm.run$S2 + allvh.eqbm.run$E2 + allvh.eqbm.run$I2) / area      # density snails of size class 2
  allvh.eqbm.run$t.3 = (allvh.eqbm.run$S3 + allvh.eqbm.run$E3 + allvh.eqbm.run$I3) / area      # density snails of size class 3
  # Estimated prevalence, using a negative binomial dist. with k = 0.2(fitted from EPLS data in nb_fit.R)
  allvh.eqbm.run$W = cov*allvh.eqbm.run$Wt + (1-cov)*allvh.eqbm.run$Wu
  allvh.eqbm.run$prev = pnbinom(2, size = 0.2, mu = allvh.eqbm.run$W, lower.tail = FALSE)   

#plot infection classes over time  
plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Snail density', ylim = c(0,max(allvh.eqbm.run$N.t)), 
     main = 'Snail Infection Classes')
  lines(allvh.eqbm.run$time, allvh.eqbm.run$I.t, col = 'red', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$E.t, col = 'orange', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$S.t, col = 'green', lwd = 2)
    legend('topright', legend = c('total', 'S', 'E', 'I'), lwd = 2, 
           col = c('black', 'green', 'orange', 'red'), cex = 0.7, bty = 'n')

#plot size classes over time
plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Snail density', ylim = c(0,max(allvh.eqbm.run$N.t)),
     main = 'Snail Size Classes')
  lines(allvh.eqbm.run$time, allvh.eqbm.run$t.1, col = 'green', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$t.2, col = 'blue', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$t.3, col = 'red', lwd = 2)
    legend('topright', legend = c('total', '1', '2', '3'), lwd = 2, col = c('black', 'green', 'blue', 'red'), cex = 0.7)
  

#plot worm burden over time  
plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$W, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Worm burden', ylim = c(0,max (allvh.eqbm.run$W)),
     main = 'Worm Burden')

#plot prevalence over time
plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$prev, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Prevalence', ylim = c(0,max (allvh.eqbm.run$prev)),
     main = 'Estimated Prevalence')

#plot infection dynamics stratified by size class
opar<-par()
par(mfrow = c(3,1), mar = c(4,3.75,1,0.4)+0.1)

plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$S1/area, type = 'l', col = 'green', lwd=2, 
     xlab = '', ylab = 'Snail density / m^2', 
     ylim = c(0,max(allvh.eqbm.run$S1/area)),
     main = 'Snail Infection Classes, size class 1')
  lines(allvh.eqbm.run$time, allvh.eqbm.run$E1/area, col = 'orange', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$I1/area, col = 'red', lwd = 2)
  legend('topright', legend = c('S', 'E', 'I'), lwd = 2, col = c('green', 'orange', 'red'), cex = 0.7)
  
plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$S2/area, type = 'l', col = 'green', lwd=2,
     xlab = '', ylab = 'Snail density / m^2', 
     ylim = c(0,max(allvh.eqbm.run$S2/area)),
     main = 'Snail Infection Classes, size class 2')
  lines(allvh.eqbm.run$time, allvh.eqbm.run$E2/area, col = 'orange', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$I2/area, col = 'red', lwd = 2)
  
plot(x = allvh.eqbm.run$time, y = allvh.eqbm.run$S3/area, type = 'l', col = 'green', lwd=2, 
     xlab = 'Time (days)', ylab = 'Snail density / m^2', 
     ylim = c(0,max(allvh.eqbm.run$S3/area)),
     main = 'Snail Infection Classes, size class 3')
  lines(allvh.eqbm.run$time, allvh.eqbm.run$E3/area, col = 'orange', lwd = 2)
  lines(allvh.eqbm.run$time, allvh.eqbm.run$I3/area, col = 'red', lwd = 2)

  par(opar) #Return to default plot settings
  