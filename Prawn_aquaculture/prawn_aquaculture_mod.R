## Modeling prawn growth and mortality in aquaculture

require(deSolve)
source("Prawn_aquaculture/macrobrachium_aquaculture_data.R")

# Model ###########
prawn_biomass = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
    
    Bp = a.p*L^b.p # Mean prawn mass, converting from length (mm) to weight (g)
    
    Bm = P*Bp # Total prawn biomass
    
    dLdt = k/(1+gam*Bm)*(linf - L) # Mean prawn length, using von Bertalanffy growth limited by density-dependent parameter gamma
    
    dPdt = -P*(muP*Bp^d + om*Bm) # Prawn abundance, subject to size- and density-dependent mortality
    
    return(list(c(dPdt,dLdt)))
  }) 
} 

# Set initial values and parameters #######
# Stocking conditions: L = 33mm ~ W = 0.5g
nstart.p = c(P = 5000, L = 33)
t.p = seq(0, 365*2, 1)
area = 1000

par.aqua=c(
  a.p = 1.21e-6,   # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii males)
  b.p = 3.43,         # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, males)
  gam = 3.5e-6,           # Density-dependent growth parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  muP = 2.21/365,       # Natural prawn mortality rate (M. volenhovenii males) from Nwosu & Wolfi 2006
  d = -0.382,           # Exponential relationship of weight with mortality, from Lorenzen 1996 (pond aquaculture)
  om = 5e-9,          # Density-dependent mortality parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010 in fit_dens_dep_params.R
  k = 1.24/365,         # Growth rate (mm/day), from Nwosu & Wolfi 2006 (M. vollenhovenii males); alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.01236667 (0.371/month)
  linf = 213.63            # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii males)
)

# Economic parameters (price estimates from Dasgupta and Tidwell)
p = 12                                           # Weighted average market price of prawns, in dollars/kg 
cost = 0.05                                          # Cost of juveniles, in dollars per
delta = -log(1-0.1)/365                          # Discount rate, equivalent to 10%/year

#Alternate growth parameter for M. rosenbergii
  k.ros <- 0.313/30
  k.vol <- 1.24/365

#run to get reference for length to weight relationship 
nstart.ref = c(P = 1000, L = 1)
p.reference = as.data.frame(ode(nstart.ref,t.p,prawn_biomass,par.aqua))
  p.reference$B = par.aqua['a.p']*p.reference$L^par.aqua['b.p']                # Mean prawn biomass, transformed from length
  p.reference$Bt = p.reference$B*p.reference$P                                                   # Total prawn biomass
 
  plot(p.reference$L, p.reference$B, type = 'l', xlab = "length (mm)", ylab = "mass (g)")
  