## Playing with scaling of area-dependent parameters in basic homogeneous predation model

require(deSolve)

area_scaling_model = function(t, n, parameters) { 
  with(as.list(parameters), {
    
    S = n[1]
    
    P.bm = (ap*(L/10)^bp)/10 # Prawn length (mm) to weight (g) conversion
    
    S.bm = as*0.8^bs # Snail length-to-weight conversion, assuming mean size = 8mm

    rps = P.bm / S.bm # Prawn-to-snail mass ratio
    
    # Attack rate as a function of biomass ratio, fit from Sokolow et al. 2014 (per m^2)
    alpha = 0.037192*rps
    
    # Adjusted attack rate, accounting for area of interest and limiting factors in the wild
    alpha_star = alpha*s/A
    
    # Handling time as a function of biomass ratio, fit from Sokolow et al. 2014
    handle = 1/(0.40450*rps)

    # Functional response
    psi = (P*alpha_star) / (1+alpha_star*handle*S)
    
    # Snail dynamics
    dSdt = f*(1-S/(Kn*A))*(S) - muN*S - psi*S
    
    
    return(list(c(dSdt)))
  }) 
}

# Set initial values and parameters
nstart = c(S = 8000)
time = seq(0, 365*2, 1)

parameters=c(
  # Prawn model parameters
  P = 50,        # Number of prawns, currently held constant
  L = 95,        # Length of prawns (mm)
  ap = 0.096868, # Allometric coefficient based on Lalrisanga et al. 2012 (for M. rosenbergii)
  bp = 3.2944,   # Allometric coefficient based on Lalrisanga et al. 2012 (for M. rosenbergii)
  
  # Snail model parameters
  as = 0.187178454, # Allometric coefficient fitted to Sanna's data on B. glabrata
  bs = 2.536764792, # Allometric coefficient fitted to Sanna's data on B. glabrata
  f = 0.16,         # Birth rate of adult snails, from Sokolow et al. 2015
  muN = 1/50,       # Natural mortality rate of snails (deaths/snail/day; assume mean lifespan = 50 days)
  Kn = 50,          # Carrying capacity of snails (per m^2, Sokolow et al. 2015)
  A = 200,          # Area of site of interest, m^2
  s = 1/10          # Scale factor limiting prawn attack rate in the wild (vs. lab conditions)
)

# Run & plot
output = as.data.frame(ode(nstart,time,area_scaling_model,parameters))
plot(output$time, output$S, type = 'l', col = 'blue', lwd=2, xlab = "Time (days)", ylab = "Snail pop.",
     ylim = c(0, 10000), main = paste('P = ', parameters['P'], ',  L = ', parameters['L'], sep = ''))

# Plot equilibrium snail abundance vs. prawn density
x = c(0:100)
s.pop = numeric(length(x))
for (i in x) {
  parameters["P"] = i
  output = as.data.frame(ode(nstart,time,area_scaling_model,parameters))
  s.pop[i+1] = tail(output$S, n = 1)
}
plot(x, s.pop, type = 'l', col = 'blue', lwd=2, xlab = "Number of prawns", ylab = "Equilibrium snail pop.",
     ylim = c(0, 10000), main = paste('L = ', parameters['L'], sep = ''))



