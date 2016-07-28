## Modeling dynamics of prawn predation on snails; no disease

## To-do list:
  # Integrate prawn growth dynamics
  # Note: attack rate scale factor changed from 1/10 in homogeneous model to 1/50 in size class model 
  #   to match expected elimination threshold; may require further exploration/tweaking

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
    alpha1 = ifelse(-log(3) + 1.1906*log(rps1) > 0, -log(3) + 1.1906*log(rps1), 0)
    alpha2 = ifelse(-log(3) + 1.1906*log(rps2) > 0, -log(3) + 1.1906*log(rps2), 0)
    alpha3 = ifelse(-log(3) + 1.1906*log(rps3) > 0, -log(3) + 1.1906*log(rps3), 0)
    
    # Adjusted attack rates, accounting for area of interest and limiting factors in the wild
    alpha_star1 = alpha1*s/A
    alpha_star2 = alpha2*s/A
    alpha_star3 = alpha3*s/A
    
    # Handling times as a function of biomass ratio, fit from Sokolow et al. 2014
    handle1 = 1/(0.38561*rps1)
    handle2 = 1/(0.38561*rps2)
    handle3 = 1/(0.38561*rps3)
    
    # Prawn dynamics, ignoring for now
    # P.bmt = P*P.bm # Total prawn biomass
    # dPdt = -P*(muP*L^-0.25 + P.bmt/phi) # Number of prawns, subject to baseline and density-dependent mortality
    # dLdt = k/(1+gam*P.bmt)*(linf - L) # Mean prawn length, with growth rate k, max length linf, crowding parameter gam
    
    # Functional response by snail size class
    psi1 = (P*alpha_star1) / (1+sum(alpha_star1*handle1*S1, alpha_star2*handle2*S2, alpha_star3*handle3*S3))
    psi2 = (P*alpha_star2) / (1+sum(alpha_star1*handle1*S1, alpha_star2*handle2*S2, alpha_star3*handle3*S3))
    psi3 = (P*alpha_star3) / (1+sum(alpha_star1*handle1*S1, alpha_star2*handle2*S2, alpha_star3*handle3*S3))
    
    # Total snail population  
    N = S1+S2+S3
    
    # Snail dynamics
    dS1dt = f*(1-N/(Kn*A))*(S2+S3) - muN1*S1 - psi1*S1 - g1*S1
    
    dS2dt = g1*S1 - muN2*S2 - psi2*S2 - g2*S2
    
    dS3dt = g2*S2 - muN3*S3 - psi3*S3
    
    
    return(list(c(dS1dt, dS2dt, dS3dt)))
  }) 
} 

# Set initial values and parameters
nstart = c(S1 = 4000, S2 = 2000, S3 = 2000)
time = seq(0, 365*2, 1)

parameters=c(
  # Prawn model parameters
  P = 50,        # Number of prawns, currently held constant
  L = 95,        # Length of prawns (mm)
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
  Kn = 50,          # Carrying capacity of snails (per m^2, Sokolow et al. 2015)
  A = 200,          # Area of site of interest, m^2
  s = 1/30          # Scale factor limiting prawn attack rate in the wild (vs. lab conditions);
                    #   decreased by factor of 3 in size class model to match expected elimination threshold
)


# Run & plot
output_ss = as.data.frame(ode(nstart,time,snail_size_predation,parameters))
plot(output_ss$time, output_ss$S1, type = 'l', col = 'red', lwd=2, xlab = "Time (days)", ylab = "Snail pop.",
     ylim = c(0, 10000), main = paste('P = ', parameters['P'], ',  L = ', parameters['L'], sep = ''))
lines(output_ss$time, output_ss$S2, col = 'green', lwd=2)
lines(output_ss$time, output_ss$S3, col = 'blue', lwd=2)
lines(output_ss$time, output_ss$S1+output_ss$S2+output_ss$S3, col = 'black', lwd=2)
legend("topright", c("S","M","L","T"), col = c("red","green","blue","black"), lty = 1, lwd = 2)

# Plot equilibrium snail abundance vs. prawn density
x = c(0:100)
s.pop = numeric(length(x))
for (i in x) {
  parameters["P"] = i
  output_ss = as.data.frame(ode(nstart,time,snail_size_predation,parameters))
  s.pop[i+1] = tail(output_ss$S1, n = 1)+tail(output_ss$S2, n = 1)+tail(output_ss$S3, n = 1)
}
plot(x, s.pop, type = 'l', col = 'blue', lwd=2, xlab = "Number of prawns", ylab = "Equilibrium snail pop.",
     ylim = c(0, 10000), main = paste('L = ', parameters['L'], sep = ''))


