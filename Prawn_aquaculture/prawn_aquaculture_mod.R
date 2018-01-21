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
    
    dPdt = -P*(muP*Bp^d + phi*Bm) # Prawn abundance, subject to size- and density-dependent mortality
    
    return(list(c(dPdt,dLdt)))
  }) 
} 

# Set initial values and parameters #######
# Stocking conditions: L = 25mm ~ W = 0.2g
nstart.p = c(P = 5000, L = 25)
t.p = seq(0, 365*2, 1)

par.aqua=c(
  a.p = 0.096868,        # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  b.p = 3.2944,          # Allometric parameter for prawn length-weight relationship, from Lalrinsanga et al. 2012 (M. rosenbergii, growout phase)
  gam = 1e-5,           # Density-dependent growth parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  muP = 0.00610958904,  # Prawn mortality at unit weight, from Lorenzen 1996 (pond aquaculture)
  d = -0.382,           # Exponential relationship of weight with mortality, from Lorenzen 1996 (pond aquaculture)
  phi = 5e-8,           # Density-dependent mortality parameter (based on biomass per hectare); informally adjusted based on Ranjeet & Kurup 2010
  k = 0.00339726,       # Growth rate (mm/day), from Nwosu & Wolfi 2006 (M. vollenhovenii); alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333
  linf = 206            # Max length (mm), from Nwosu & Wolfi 2006 (M. vollenhovenii)
)

# Economic parameters (price estimates from Tamil Nadu Agricultural University, http://agritech.tnau.ac.in/fishery/fish_freshwaterprawn.html)
p = 12                                           # Weighted average market price of prawns, in rupees/kg 
c = 0.1                                           # Cost of post-larvae, in rupees/1000 PL
delta = -log(1-0.1)/365                           # Discount rate, equivalent to 10%/year


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
