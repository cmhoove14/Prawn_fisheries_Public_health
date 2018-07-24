## Modeling prawn growth and mortality in aquaculture

require(deSolve)

# Model ###########
prawn_biomass = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
    
    Bp = (a.p/10)*(L/10)^b.p # Mean prawn mass, converting from length (mm) to weight (g)
    
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

par.aqua=c(
  a.p = 0.096868,        # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  b.p = 3.2944,          # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  gam = 6e-6,           # Density-dependent growth parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  muP = 0.00610958904,  # Prawn mortality at unit weight, from Lorenzen 1996 (pond aquaculture); informally adjusted based on Ranjeet & Kurup 2010 in fit_dens_dep_params.R
  d = -0.382,           # Exponential relationship of weight with mortality, from Lorenzen 1996 (pond aquaculture)
  om = 5e-9,            # Density-dependent mortality parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010 in fit_dens_dep_params.R
  k = 0.00339726,       # Growth rate (mm/day), from Nwosu & Wolfi 2006 (M. vollenhovenii); alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333
  linf = 206            # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii)
)

# Economic parameters (price estimates from Dasgupta and Tidwell)
p = 12                                           # Weighted average market price of prawns, in dollars/kg 
c = 0.1                                          # Cost of post-larvae, in dollars per
delta = -log(1-0.1)/365                          # Discount rate, equivalent to 10%/year


#run to get optimization parameters
p.run = as.data.frame(ode(nstart.p,t.p,prawn_biomass,par.aqua))
  p.run$B = ((par.aqua['a.p']*(p.run$L/10)^par.aqua['b.p'])/10)                # Mean prawn biomass, transformed from length
  p.run$Bt = p.run$B*p.run$P                                                   # Total prawn biomass

  start.mass.kg = p.run$Bt[p.run$time==0]/1000
  harvest.mass.kg = max(p.run$Bt)/1000
  harvest.b = p.run$B[p.run$Bt==max(p.run$Bt)] 
  harvest.p = p.run$P[p.run$Bt==max(p.run$Bt)] 
  harvest.l = p.run$L[p.run$Bt==max(p.run$Bt)] 
  harvest.time = p.run$time[p.run$Bt==max(p.run$Bt)]

#run to get reference for length to weight relationship 
nstart.ref = c(P = 10000, L = 1)
p.reference = as.data.frame(ode(nstart.ref,t.p,prawn_biomass,par.aqua))
  p.reference$B = ((par.aqua['a.p']*(p.reference$L/10)^par.aqua['b.p'])/10)                # Mean prawn biomass, transformed from length
  p.reference$Bt = p.reference$B*p.reference$P                                                   # Total prawn biomass
 