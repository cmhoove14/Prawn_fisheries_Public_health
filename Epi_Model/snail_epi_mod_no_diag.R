## Epidemiological model including snail size classes; no predation

## To-do list:
  # Continue tuning infection parameters, draw estimates from literature as necessary
  # May want to consider a few aspects of snail dynamics further, e.g.: 
  # Should exposed snails be able to reproduce at a reduced rate?
  # Should growth rates differ by infection status?
  # Think a little more about diagonal transitions
  # Reconsider area scaling
  # Come up with R0 expression?

require(deSolve)
require(tidyverse)

#mating function
fx<-function(x, mean.worm, clump){
  alpha<-(mean.worm)/(clump+mean.worm)
  (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
}
phi_Wk<-function(W, phi){
  alpha<-W/(W+phi)
  1-( (1-alpha)^(phi+1) * (integrate(fx, 0, 2*pi, W, phi, stop.on.error = F)$value) /(2*pi)  )
}

# Set initial values and parameters #######
area = 1000 #1 square kiolmeter = 1000m^2
nstart.sn = c(S1 = 5*area, S2 = 0, S3 = 0, 
              E1 = 0, E2 = 0, E3 = 0, 
              I1 = 0, I2 = 0, I3 = 0, Wt = 10, Wu = 10)
t.sn = seq(0,365*300,10)
cov = 0.8

par.snails=c(
  ## Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5 * area,          # Human population at site of interest
  
  ## Reproductive parameters
  f = 0.26,          # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection; more like a recruitment rate)
  Kn = 50*area,      # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.5,           # Fraction of exposed snails that can reproduce, from Sokolow et al. 2015
  
  ## Snail mortality parameters
  muN1 = 1/50,       # Natural mortality rate of small snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN2 = 1/75,       # Natural mortality rate of medium snails (deaths/snail/day; assume mean lifespan = 75 days)
  muN3 = 1/100,      # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 100 days)
  muI = 1/10,        # Mortality rate of infected snails

  ## Predation parameters
  psi1 = 0,          # Predation rate of prawns on small snails, ignored in this model
  psi2 = 0,          # Predation rate of prawns on medium snails, ignored in this model
  psi3 = 0,          # Predation rate of prawns on large snails, ignored in this model
  
  ## Snail growth parameters
  g1 = 1/37,         # Growth rate of small snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  g2 = 1/62,         # Growth rate of medium snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  
  ## Infection parameters
  beta = 4e-6,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day); from Sokolow et al. 2015 
  m = 0.8,           # Miracidial shedding rate per adult female worm divided by miracidial mortality; from Sokolow et al. 2015
  sigma = 1/50,      # Latent period for exposed snails (infectious snails/exposed snail/day); from Sokolow et al. 2015 
  lambda = 0.005,  # Snail-to-human infection probability scaled to 1 m^2 (composite including cercarial shedding, mortality, infection, survival to patency) from Sokolow et al. 2015
  theta1 = 167.8/127.8,         # Scale factor describing increase in cercarial shedding rate in larger (size class 2) snails; from Chu & Dawood 1970 (estimated to be between 2 and 10)
  theta2 = 1006.8/127.8,         # Scale factor describing increase in cercarial shedding rate in larger (size class 2) snails; from Chu & Dawood 1970 (estimated to be between 2 and 10)

  phi = 0.2,           # clumping parameter of the negative binomial distribution used in the mating probability function
  
  ## Schisto mortality parameters
  muW = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  muH = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)




#Model with no diagonal transitions, all horizontal and vertical transitions #####
snail_epi_allvh = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S1=n[1]
    S2=n[2]
    S3=n[3]
    E1=n[4]
    E2=n[5]
    E3=n[6]
    I1=n[7]
    I2=n[8]
    I3=n[9]
    Wt=n[10]
    Wu=n[11]
    
    W = cov*Wt + (1-cov)*Wu
    
    mate = phi_Wk(W = W, phi = phi)  #Mating probability
    S = S1+S2+S3    # Total susceptible snails
    E = E1+E2+E3    # Total exposed snails
    I = I1+I2+I3       # Total infected snails
    N = S+E+I       # Total number of snails
    M = m*0.5*W*mate # Miracidial density per person as a function of mean worm burden (W), miracidial shedding rate (m), and mating probability
    # rho = g2/sigma  # Fraction of E2 snails that transition directly to I3
    
    dS1dt = f*(1-N/Kn)*(S2+S3+z*(E2+E3)) - muN1*S1 - psi1*S1 - g1*S1 - (beta)*M*H*S1
    
    dS2dt = g1*S1 - muN2*S2 - psi2*S2 - g2*S2 - (beta)*M*H*S2
    
    dS3dt = g2*S2 - muN3*S3 - psi3*S3 - (beta)*M*H*S3
    
    dE1dt = (beta)*M*H*S1 - muN1*E1 - psi1*E1 - sigma*E1 - g1*E1
    
    dE2dt = (beta)*M*H*S2 + g1*E1 - muN2*E2 - psi2*E2 - sigma*E2 - g2*E2
    
    dE3dt = (beta)*M*H*S3 + g2*E2 - muN3*E3 - psi3*E3 - sigma*E3
    
    dI1dt = sigma*E1 - (muN1+muI)*I1 - psi2*I1 - g1*I1

    dI2dt = sigma*E2 + g1*I1 - (muN2+muI)*I2 - psi2*I2 - g2*I2
    
    dI3dt = sigma*E3 + g2*I2 - (muN3+muI)*I3 - psi3*I3
    
    dWtdt = lambda*I1/A + theta1*lambda*I2/A + theta2*lambda*I3/A - (muW + muH)*Wt
    
    dWudt = lambda*I1/A + theta1*lambda*I2/A + theta2*lambda*I3/A - (muW + muH)*Wu
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI1dt, dI2dt, dI3dt, dWtdt, dWudt)))
  }) 
} 

#Run to equilibrium to eqbm values of state parameters
allvh.eqbm.run = as.data.frame(ode(nstart.sn,t.sn,snail_epi_allvh,par.snails))
  allvh.eqbm = allvh.eqbm.run[dim(allvh.eqbm.run)[1],]
