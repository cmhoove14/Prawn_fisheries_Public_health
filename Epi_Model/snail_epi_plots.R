## Epidemiological model including snail size classes; no predation

## To-do list:
  # Continue tuning infection parameters, draw estimates from literature as necessary
  # May want to consider a few aspects of snail dynamics further, e.g.: 
      # Should exposed snails be able to reproduce at a reduced rate?
      # Should growth rates differ by infection status?
      # Think a little more about diagonal transitions
      # Reconsider area scaling
  # Come up with R0 expression?

source('Epi_Model/snail_epi_mod.R')

# Run & plot #########
output_s.epi = as.data.frame(ode(nstart.sn,t.sn,snail_epi,par.snails))

  snailepiq<-output_s.epi[dim(output_s.epi)[1],] #save equilibrium values

  output_s.epi$S.t = (output_s.epi$S1 + output_s.epi$S2 + output_s.epi$S3) / area        # density susceptible snails
  output_s.epi$E.t = (output_s.epi$E1 + output_s.epi$E2 + output_s.epi$E3) / area        # density exposed snails 
  output_s.epi$I.t = (output_s.epi$I2 + output_s.epi$I3 ) / area                         # density infected snails
  output_s.epi$N.t = (output_s.epi$S.t + output_s.epi$E.t + output_s.epi$I.t)            # density snails
  output_s.epi$t.1 = (output_s.epi$S1 + output_s.epi$E1) / area                          # density snails of size class 1
  output_s.epi$t.2 = (output_s.epi$S2 + output_s.epi$E2 + output_s.epi$I2) / area      # density snails of size class 2
  output_s.epi$t.3 = (output_s.epi$S3 + output_s.epi$E3 + output_s.epi$I3) / area      # density snails of size class 3
  # Estimated prevalence, using a negative binomial dist. with k = 0.2(fitted from EPLS data in nb_fit.R)
  output_s.epi$prev = pnbinom(2, size = 0.2, mu = output_s.epi$W, lower.tail = FALSE)   

#plot infection classes over time  
plot(x = output_s.epi$time, y = output_s.epi$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Snail density', ylim = c(0,max(output_s.epi$N.t)), 
     main = 'Snail Infection Classes')
  lines(output_s.epi$time, output_s.epi$I.t, col = 'red', lwd = 2)
  lines(output_s.epi$time, output_s.epi$E.t, col = 'orange', lwd = 2)
  lines(output_s.epi$time, output_s.epi$S.t, col = 'green', lwd = 2)
    legend('topright', legend = c('total', 'S', 'E', 'I'), lwd = 2, 
           col = c('black', 'green', 'orange', 'red'), cex = 0.7, bty = 'n')

#plot size classes over time
plot(x = output_s.epi$time, y = output_s.epi$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Snail density', ylim = c(0,max(output_s.epi$N.t)),
     main = 'Snail Size Classes')
  lines(output_s.epi$time, output_s.epi$t.1, col = 'green', lwd = 2)
  lines(output_s.epi$time, output_s.epi$t.2, col = 'blue', lwd = 2)
  lines(output_s.epi$time, output_s.epi$t.3, col = 'red', lwd = 2)
    legend('topright', legend = c('total', '1', '2', '3'), lwd = 2, col = c('black', 'green', 'blue', 'red'), cex = 0.7)
  

#plot worm burden over time  
plot(x = output_s.epi$time, y = output_s.epi$W, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Worm burden', ylim = c(0,max (output_s.epi$W)),
     main = 'Worm Burden')

#plot prevalence over time
plot(x = output_s.epi$time, y = output_s.epi$prev, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Prevalence', ylim = c(0,max (output_s.epi$prev)),
     main = 'Estimated Prevalence')

#plot infection dynamics stratified by size class
opar<-par()
par(mfrow = c(3,1), mar = c(4,3.75,1,0.4)+0.1)

plot(x = output_s.epi$time, y = output_s.epi$S1/area, type = 'l', col = 'green', lwd=2, 
     xlab = '', ylab = 'Snail density / m^2', 
     ylim = c(0,max(output_s.epi$S1/area)),
     main = 'Snail Infection Classes, size class 1')
  lines(output_s.epi$time, output_s.epi$E1/area, col = 'orange', lwd = 2)
  legend('topright', legend = c('S', 'E', 'I'), lwd = 2, col = c('green', 'orange', 'red'), cex = 0.7)
  
plot(x = output_s.epi$time, y = output_s.epi$S2/area, type = 'l', col = 'green', lwd=2,
     xlab = '', ylab = 'Snail density / m^2', 
     ylim = c(0,max(output_s.epi$S2/area)),
     main = 'Snail Infection Classes, size class 2')
  lines(output_s.epi$time, output_s.epi$E2/area, col = 'orange', lwd = 2)
  lines(output_s.epi$time, output_s.epi$I2/area, col = 'red', lwd = 2)
  
plot(x = output_s.epi$time, y = output_s.epi$S3/area, type = 'l', col = 'green', lwd=2, 
     xlab = 'Time (days)', ylab = 'Snail density / m^2', 
     ylim = c(0,max(output_s.epi$S3/area)),
     main = 'Snail Infection Classes, size class 3')
  lines(output_s.epi$time, output_s.epi$E3/area, col = 'orange', lwd = 2)
  lines(output_s.epi$time, output_s.epi$I3/area, col = 'red', lwd = 2)

  par(opar) #Return to default plot settings