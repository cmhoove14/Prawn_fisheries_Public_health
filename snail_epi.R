## Epidemiological model including snail size classes; no predation

## To-do list:
  # Draw parameter estimates from literature
  # Make sure snail dynamics between size and infection classes are appropriate, e.g.: 
      # Should small snails be able to get infected?
      # Should exposed snails be able to reproduce at a reduced rate?
      # Think a little more about diagonal transitions
  # Come up with R0 expression

require(deSolve)

snail_epi = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S1=n[1]
    S2=n[2]
    S3=n[3]
    E1=n[4]
    E2=n[5]
    E3=n[6]
    I2=n[7]
    I3=n[8]
    C=n[9]
    W=n[10]
    
    S = S1+S2+S3  # Total susceptible snails
    E = E1+E2+E3  # Total exposed snails
    I = I2+I3     # Total infected snails
    N = S+E+I     # Total number of snails
    M = m*0.5*W   # Miracidial density as a function of mean worm burden (W) and miracidial shedding rate (m)
                  # Note: add mating function?
    
    dS1dt = f*(1-N/(Kn*A))*(S2+S3+z*(E2+E3)) - muN1*S1 - psi1*S1 - g1*S1 - (beta/A)*M*H*S1
    
    dS2dt = g1*S1 - muN2*S2 - psi2*S2 - g2*S2 - (beta/A)*M*H*S2
    
    dS3dt = g2*S2 - muN3*S3 - psi3*S3 - (beta/A)*M*H*S3
    
    dE1dt = (beta/A)*M*H*S1 - muN1*E1 - psi1*E1 - g1*E1 - sigma*X*E1
    
    dE2dt = (beta/A)*M*H*S2 + g1*E1 - muN2*E2 - psi2*E2 - g2*E2 - sigma*E2
    
    dE3dt = (beta/A)*M*H*S3 + g2*E2 - muN3*E3 - psi3*E3 - sigma*E3
    
    dI2dt = sigma*X*E1 + sigma*(1-rho)*E2 - (muN2 + muI)*I2 - psi2*I2 - g2*I2
    
    dI3dt = sigma*rho*E2 + sigma*E3 + g2*I2 - (muN3 + muI)*I3 - psi3*I3
    
    dCdt = theta2*I2 + theta3*I3 - muC*C
    
    dWdt = lambda*C - (muW + muH)*W
    
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI2dt, dI3dt, dCdt, dWdt)))
  }) 
} 

# Set initial values and parameters
nstart = c(S1 = 1, S2 = 0, S3 = 0, E1 = 0, E2 = 0, E3 = 0, I2 = 0, I3 = 0, C = 0, W = 2)
time = seq(0,365*2,1)

parameters=c(
  ## Location parameters
  A = 200,           # Area of site of interest, m^2
  H = 1000,          # Human population at site of interest
  
  ## Reproductive parameters
  f = 0.26,          # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection; more like a recruitment rate)
  Kn = 50,           # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.5,           # Fraction of exposed snails that can reproduce, from Sokolow et al. 2015
  
  ## Snail mortality parameters
  muN1 = 1/40,       # Natural mortality rate of small snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN2 = 1/50,       # Natural mortality rate of medium snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN3 = 1/60,       # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 50 days)
  muI = 1/10,        # Additional mortality rate of shedding snails as a result of infection, from Sokolow et al. 2015
  
  ## Predation parameters
  psi1 = 0,          # Predation rate of prawns on small snails, ignored in this model
  psi2 = 0,          # Predation rate of prawns on medium snails, ignored in this model
  psi3 = 0,          # Predation rate of prawns on large snails, ignored in this model
  
  ## Snail growth parameters
  g1 = 1/37,         # Growth rate of small snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  g2 = 1/62,         # Growth rate of medium snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  
  ## Infection parameters
  beta = 8e-4,       # Human-to-snail infection probability (infected snails/miracidia/snail/day), from Sokolow et al. 2015, scaled to 1 m^2
  m = 0.8,           # Miracidial shedding rate per adult female worm divided by miracidial mortality, from Sokolow et al. 2015
  sigma = 1/50,      # Latent period for exposed snails (infectious snails/exposed snail/day), from Sokolow et al. 2015
  X = 0,             # Fraction of exposed small snails that convert directly to medium shedding; no source
  rho = 0,           # Fraction of exposed medium snails that convert directly to large shedding; no source
  theta2 = 500,      # Cercarial shedding rate for medium infected snails (cercariae/infectious medium snail/day); no source
  theta3 = 1000,     # Cercarial shedding rate for large infected snails (cercariae/infected large snail/day); no source
  lambda = 1e-4,     # Probability of cercaria successfully infecting human and maturing; no source
  
  ## Schisto mortality parameters
  muC = 0.99,        # Cercarial mortality rate (very small chance of cercariae lasting >1 day); no source
  muW = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  muH = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)


# Run & plot
output_s.epi=as.data.frame(ode(nstart,time,snail_epi,parameters))
output_s.epi$S.t = output_s.epi$S1 + output_s.epi$S2 + output_s.epi$S3    # Total susceptible snails
output_s.epi$E.t = output_s.epi$E1 + output_s.epi$E2 + output_s.epi$E3    # Total exposed snails 
output_s.epi$I.t = output_s.epi$I2 + output_s.epi$I3                      # Total infected snails
output_s.epi$N.t = output_s.epi$S.t + output_s.epi$E.t + output_s.epi$I.t # Total snails
output_s.epi$t.1 = output_s.epi$S1 + output_s.epi$E1                      # Total snails of size class 1
output_s.epi$t.2 = output_s.epi$S2 + output_s.epi$E2 + + output_s.epi$I2  # Total snails of size class 2
output_s.epi$t.3 = output_s.epi$S3+ output_s.epi$E3 + + output_s.epi$I3   # Total snails of size class 3
  
plot(x = output_s.epi$time, y = output_s.epi$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output_s.epi$N.t)),
     main = 'Snail Size Classes')
  lines(output_s.epi$time, output_s.epi$t.1, col = 'green', lwd = 2)
  lines(output_s.epi$time, output_s.epi$t.2, col = 'blue', lwd = 2)
  lines(output_s.epi$time, output_s.epi$t.3, col = 'red', lwd = 2)
  legend('topright', legend = c('total', '1', '2', '3'), lwd = 2, col = c('black', 'green', 'blue', 'red'), cex = 0.7)
  
plot(x = output_s.epi$time, y = output_s.epi$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output_s.epi$N.t)), 
     main = 'Snail Infection Classes')
  lines(output_s.epi$time, output_s.epi$I.t, col = 'red', lwd = 2)
  lines(output_s.epi$time, output_s.epi$E.t, col = 'orange', lwd = 2)
  lines(output_s.epi$time, output_s.epi$S.t, col = 'green', lwd = 2)
  legend('topright', legend = c('total', 'S', 'E', 'I'), lwd = 2, col = c('black', 'green', 'orange', 'red'), cex = 0.7)
  
plot(x = output_s.epi$time, y = output_s.epi$C, type = 'l', col = 'orange', lwd=2, xlab = 'Time (days)', 
     ylab = 'Cercariae', ylim = c(0,max (output_s.epi$C)),
     main = 'Cercarial Density')

plot(x = output_s.epi$time, y = output_s.epi$W, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Worm burden', ylim = c(0,max (output_s.epi$W)),
     main = 'Worm Burden')

