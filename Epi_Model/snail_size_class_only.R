## Modeling base snail population dynamics to set growth and mortality parameters; 
## no predation or disease

## To-do list:
  # Size-dependent mortality in snail classes? (simple approximation looks decent)

require(deSolve)

snail_size <- function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S1 = n[1]
    S2 = n[2]
    S3 = n[3]
  
    N = S1 + S2 + S3
    
    dS1dt = f*(1-N/Kn)*(S2+S3) - muN1*S1 - psi1*S1 - g1*S1
    
    dS2dt = g1*S1 - muN2*S2 - psi2*S2 - g2*S2
    
    dS3dt = g2*S2 - muN3*S3 - psi3*S3
    
    
    return(list(c(dS1dt, dS2dt, dS3dt)))
  }) 
} 

# Set initial values and parameters
nstart = c(S1 = 0, S2 = 0.1, S3 = 0)
time = seq(0, 365*2, 1)

# List snail model parameters
parameters = c(
    f = 0.26,    # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection; more like a recruitment rate)
    Kn = 50,     # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
    muN1 = 1/40, # Natural mortality rate of small snails (deaths/snail/day; assume mean lifespan = 50 days)
    muN2 = 1/50, # Natural mortality rate of medium snails (deaths/snail/day; assume mean lifespan = 50 days)
    muN3 = 1/60, # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 50 days)
    g1 = 1/37,   # Growth rate of small snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C) 
    g2 = 1/62,   # Growth rate of medium snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C) 
    psi1 = 0,    # TBD: predation rate of prawn cohort on smallest snails
    psi2 = 0,    # TBD: predation rate of prawn cohort on medium snails
    psi3 = 0     # TBD: predation rate of prawn cohort on largest snails
)


# Run & plot
output_ss = as.data.frame(ode(nstart,time,snail_size,parameters))
output_ss$S.t = output_ss$S1 + output_ss$S2 + output_ss$S3 # Calculate total snail population density
eqbm_ss = output_ss[max(time),]
  
par(mfrow = c(1,1))  
plot(output_ss$time, output_ss$S.t, type = 'l', col = 'black', lwd=2, ylim = c(-2, 45), #xlim = c(0,305),
     ylab = 'state variables', xlab = 'time',
     main = paste('Kn=', parameters['Kn'], '  muN1=', parameters['muN1'], '  muN2=', parameters['muN2'],
                  '\n  muN3=', parameters['muN3'], '  f=', parameters['f'], sep = ''))
  lines(output_ss$time, output_ss$S1, type = 'l', col = 'green', lwd=2)
  lines(output_ss$time, output_ss$S2, type = 'l', col = 'blue', lwd=2)
  lines(output_ss$time, output_ss$S3, type = 'l', col = 'red', lwd=2)
  abline(v = 365, lty = 2)
  legend('topright', legend = c('S1', 'S2', 'S3', 'S.t'), lwd = 2, col = c('green', 'blue', 'red', 'black'))
    
# Simulate and plot single equilibrium-level cohort
nstart = c(S1 = eqbm_ss$S.t, S2 = 0, S3 = 0)
time = seq(0,365*2,1)
    
parameters['f'] = 0 # Stop recruitment
    
output_ss = as.data.frame(ode(nstart,time,snail_size,parameters))
output_ss$S.t = output_ss$S1 + output_ss$S2 + output_ss$S3 # Calculate total snail population density
eqbm_ss = output_ss[max(time),]
    
par(mfrow = c(1,1))  
plot(output_ss$time, output_ss$S.t, type = 'l', col = 'black', lwd=2, ylim = c(-2, 45), xlim = c(0,305),
     ylab = 'state variables', xlab = 'time',
     main = paste('Kn=', parameters['Kn'], '  muN1=', parameters['muN1'], '  muN2=', parameters['muN2'],
                  '\n  muN3=', parameters['muN3'], '  f=', parameters['f'], sep = ''))
  lines(output_ss$time, output_ss$S1, type = 'l', col = 'green', lwd=2)
  lines(output_ss$time, output_ss$S2, type = 'l', col = 'blue', lwd=2)
  lines(output_ss$time, output_ss$S3, type = 'l', col = 'red', lwd=2)
  abline(v = 365, lty=2)
  legend('topright', legend = c('S1', 'S2', 'S3', 'S.t'), lwd=2, col = c('green', 'blue', 'red', 'black'))
    
# Calculate PDF f(t) = -dS(t)/dt (i.e. numerically approximate derivative of survival curve)
vect = output_ss$S.t
deriv = numeric(length(vect)-1)
for (i in 1:(length(vect)-1)) {
  deriv[i] = vect[i] - vect[i+1]
}

deriv2 = deriv/output_ss$S.t[1]
deriv2 = deriv2[-(length(deriv2))]
sum(deriv2) 

x = 0:(length(deriv2)-1)

plot(x = x, y = deriv2, type = 'l')

# Calculate life expectancy
sum(x*deriv2)   

