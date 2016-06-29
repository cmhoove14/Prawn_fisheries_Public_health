## Modeling dynamics of prawn predation on snails; no disease

## To-do list:
  # Scale attack rate, density-dependent parameters appropriately by area
  # Vary prawn size, eventually include growth dynamics
  # Look for snail elimination thresholds (plot N vs. P)

require(deSolve)

snail_size_predation = function(t, n, parameters) { 
  with(as.list(parameters), {
    
    S1 = n[1]
    S2 = n[2]
    S3 = n[3]
    
    P.bm = (ap*(L/10)^bp)/10 # Prawn length (mm) to weight (g) conversion
    
    # Snail length-to-weight conversion
    S1.bm = as*0.4^bs # 4mm size class
    S2.bm = as*0.8^bs # 8mm size class
    S3.bm = as*1.2^bs # 12mm size class
    
    rps1 = P.bm / S1.bm # Prawn-to-snail mass ratio for size class 1
    rps2 = P.bm / S2.bm # Prawn-to-snail mass ratio for size class 2
    rps3 = P.bm / S3.bm # Prawn-to-snail mass ratio for size class 3
      
    # Attack rates as a function of biomass ratio, fit from Sokolow et al. 2014
    alpha1 = 0.037192*rps1
    alpha2 = 0.037192*rps2
    alpha3 = 0.037192*rps3
    
    # Handling times as a function of biomass ratio, fit from Sokolow et al. 2014
    handle1 = 1/(0.40450*rps1)
    handle2 = 1/(0.40450*rps2)
    handle3 = 1/(0.40450*rps3)
    
    # Prawn dynamics, ignoring for now
    # P.bmt = P*P.bm # Total prawn biomass
    # dPdt = -P*(muP*L^-0.25 + P.bmt/phi) # Number of prawns, subject to baseline and density-dependent mortality
    # dLdt = k/(1+gam*P.bmt)*(linf - L) # Mean prawn length, with growth rate k, max length linf, crowding parameter gam
    
    # Functional response by snail size class
    psi1 = (P*alpha1) / (1+sum(alpha1*handle1*S1, alpha2*handle2*S2, alpha3*handle3*S3))
    psi2 = (P*alpha2) / (1+sum(alpha1*handle1*S1, alpha2*handle2*S2, alpha3*handle3*S3))
    psi3 = (P*alpha3) / (1+sum(alpha1*handle1*S1, alpha2*handle2*S2, alpha3*handle3*S3))
    
    # Total snail population  
    N = S1+S2+S3
    
    # Snail dynamics
    dS1dt = f*(1-N/Kn)*(S2+S3) - muN1*S1 - psi1*S1 - g1*S1
    
    dS2dt = g1*S1 - muN2*S2 - psi2*S2 - g2*S2
    
    dS3dt = g2*S2 - muN3*S3 - psi3*S3
    
    
    return(list(c(dS1dt, dS2dt, dS3dt)))
  }) 
} 

# Set initial values and parameters
nstart = c(S1 = 5000, S2 = 2000, S3 = 1000)
time = seq(0, 365*2, 1)

parameters=c(
  # Prawn model parameters
  P = 50,        # Number of prawns, currently held constant
  L = 75,        # Length of prawns (mm)
  ap = 0.096868, # Allometric coefficient based on Lalrisanga et al. 2012 (for M. rosenbergii)
  bp = 3.2944,   # Allometric coefficient based on Lalrisanga et al. 2012 (for M. rosenbergii)
  # gam = 1e-6, # Crowding parameter, reduces growth rate at high density
  # muP = 0.006136986, # Baseline prawn mortality rate
  # phi = 40000000, # Biomass-based density dependence parameter, per hactare
  # k = 0.00339726, # Growth rate (mm/day)
  # linf = 206, # Max length (mm)
    
  # Snail model parameters
  as = 0.187178454, # Allometric coefficient fitted to Sanna's data on B. glabrata
  bs = 2.536764792, # Allometric coefficient fitted to Sanna's data on B. glabrata
  f = 0.26,         # Birth rate of adult snails, from Sokolow et al. 2015 scaled by 1/0.6 to account for immature snails
  muN1 = 1/40,      # Natural mortality rate of small snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN2 = 1/50,      # Natural mortality rate of medium snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN3 = 1/60,      # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 50 days)
  g1 = 1/37,        # Growth rate of small snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C) 
  g2 = 1/62,        # Growth rate of medium snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C) 
  Kn = 10000        # Carrying capacity of snails (per 200m^2 site, Sokolow et al. 2015) 
)


# Run & plot
output_ss = as.data.frame(ode(nstart,time,snail_size_predation,parameters))
plot(output_ss$time, output_ss$S1, type = 'l', col = 'red', lwd=2, xlab = "Time (days)", ylab = "Snail pop.",
     ylim = c(0, 5000))
lines(output_ss$time, output_ss$S2, col = 'green', lwd=2)
lines(output_ss$time, output_ss$S3, col = 'blue', lwd=2)
legend("topright", c("S","M","L"), col = c("red","green","blue"), lty = 1, lwd = 2)



