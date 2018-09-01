
# Set initial values and parameters #######
area = 1000 #1 square kiolmeter = 1000m^2
nstart.sn = c(S1 = 5*area, S2 = 5*area, S3 = 5*area, 
              E1 = 5*area, E2 = 5*area, E3 = 5*area, 
              I1 = 5*area, I2 = 5*area, I3 = 5*area, 
              Wt = 50, Wu = 50)

t.sn = seq(0,365*300,10)
cvrg = 0.75  #MDA coverage
eff = 0.85 #MDA efficacy 

# Parameters used in the epidemiological model 
par.snails.imm=c(
  ## Location parameters
  A = area,          # Area of site of interest, m^2
  H = 1.5 * area,          # Human population at site of interest
  
  ## Reproductive parameters
  f = 0.16,          # Birth rate of adult snails (snails/reproductive snail/day, including survival to detection; more like a recruitment rate)
  Kn = 50*area,      # Carrying capacity of snails (snails/m^2), from Sokolow et al. 2015
  z = 0.5,           # Fraction of exposed snails that can reproduce, from Sokolow et al. 2015
  
  ##immigration parameters
  iota = 0,       # % of snails in each infection class immigrate from other site
  siteS1 = 0,        # Density of snails in each class introduced via immigration
  siteS2 = 0,
  siteS3 = 0,
  siteE1 = 0,
  siteE2 = 0,
  siteE3 = 0,
  siteI1 = 0,
  siteI2 = 0,
  siteI3 = 0,
  
  ## Snail mortality parameters
  muN1 = 1/50,       # Natural mortality rate of small snails (deaths/snail/day; assume mean lifespan = 50 days)
  muN2 = 1/75,       # Natural mortality rate of medium snails (deaths/snail/day; assume mean lifespan = 75 days)
  muN3 = 1/100,      # Natural mortality rate of large snails (deaths/snail/day; assume mean lifespan = 100 days)
  muI = 1/10,        # Mortality rate of infected snails

  ## Predation parameters
  psi1 = 0,          # Predation rate of prawns on small snails, ignored in this model
  psi2 = 0,          # Predation rate of prawns on medium snails, ignored in this model
  psi3 = 0,          # Predation rate of prawns on large snails, ignored in this model
  nfr = 2,           # Exponent of functional response (=1 for type 2 response, =2 for type 3 response)

  ## Snail growth parameters
  g1 = 1/37,         # Growth rate of small snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  g2 = 1/62,         # Growth rate of medium snails (size class transition rate, in terms of days to grow 4mm; adapted from McCreesh et al. 2014, assuming water temp. of 25 C)
  
  ## Infection parameters
  beta = 4e-7,       # Human-to-snail infection probability in reference area (infected snails/miracidia/snail/day); from Sokolow et al. 2015 
  m = 0.8,           # Miracidial shedding rate per adult female worm divided by miracidial mortality; from Sokolow et al. 2015
  sigma = 1/30,      # Latent period for exposed snails (infectious snails/exposed snail/day); from Sokolow et al. 2015 
  lambda = 7.5e-6,  # Snail-to-human infection probability scaled to 1 m^2 (composite including cercarial shedding, mortality, infection, survival to adult worm) from Sokolow et al. 2015
  theta1 = 167.8/127.8,         # Scale factor describing increase in cercarial shedding rate in larger (size class 2) snails; from Chu & Dawood 1970 (estimated to be between 2 and 10)
  theta2 = 1006.8/127.8,         # Scale factor describing increase in cercarial shedding rate in larger (size class 3) snails; from Chu & Dawood 1970 (estimated to be between 2 and 10)

  phi = 0.08,           # clumping parameter of the negative binomial distribution used in the mating probability function
  
  ## Schisto mortality parameters
  muW = 1/(3.3*365), # Natural mortality rate of adult worms in humans, assuming average lifespan of 3.3 years, from Sokolow et al. 2015
  muH = 1/(60*365)   # Natural mortality of humans (contributing to worm mortality), assuming average lifespan of 60 years, from Sokolow et al. 2015
)
