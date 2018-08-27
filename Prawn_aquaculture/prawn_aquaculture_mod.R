## Modeling prawn growth and mortality in aquaculture

require(deSolve)
source("Prawn_aquaculture/macrobrachium_aquaculture_data.R")

# Model ###########
prawn_biomass = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
    
    Bp = 10^a.p*(L/10)^b.p  # a.p*L^b.p # Mean prawn mass, converting from length (cm) to weight (g)
    
    Bm = P*Bp # Total prawn biomass
    
    dLdt = k/(1+gam*Bm)*(linf - L) # Mean prawn length, using von Bertalanffy growth limited by density-dependent parameter gamma
    
    dPdt = -P*(muP*Bp^d + om*Bm) # Prawn abundance, subject to size- and density-dependent mortality
    
    return(list(c(dPdt,dLdt)))
  }) 
} 

# Set initial values and parameters #######
# Stocking conditions: L = 40mm ~ W = 0.4g
nstart.p = c(P = 5000, L = 40)
t.p = seq(0, 365*2, 1)
area = 1000

par.aqua=c(
  a.p = -2.6132,   # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii males)
  b.p = 3.5502,      # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, males)
  gam = 5e-6,    # Density-dependent growth parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  muP = 2.21/365,  # Natural prawn mortality rate (M. volenhovenii males) from Nwosu & Wolfi 2006
  d = -0.382,      # Exponential relationship of weight with mortality, from Lorenzen 1996 (pond aquaculture)
  om = 5.5e-9,       # Density-dependent mortality parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010 in fit_dens_dep_params.R
  k = 3.19/365,    # Growth rate (mm/day) of M. vollenhovenii
  linf = 213.63,            # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii males)
  k.ros = 0.371/30 # Growth rate (mm/day) of M. rosenbergii
)

# Economic parameters (price estimates from Dasgupta and Tidwell)
price = 12                                       # Weighted average market price of prawns, in dollars/kg 
cost = 0.10                                      # Cost of juveniles, in dollars per
delta = -log(1-0.03)/365                         # Discount rate, equivalent to 3%/year

#Alternate growth parameters for M. rosenbergii and M. vollenhovenii
  k.ros <- 0.371/30    # From Sampaio and Valenti
  k.vol <- 3.19/365    # From Gabche and Hockey 1995

#run to get reference for length to weight relationship 
nstart.ref = c(P = 1000, L = 1)
p.reference = as.data.frame(ode(nstart.ref,t.p,prawn_biomass,par.aqua))
  p.reference$B = 10^par.aqua['a.p']*(p.reference$L/10)^par.aqua['b.p']                # Mean prawn biomass, transformed from length
  p.reference$Bt = p.reference$B*p.reference$P                                                   # Total prawn biomass
 