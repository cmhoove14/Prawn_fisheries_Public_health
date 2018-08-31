# Parameters to use in the prawn aquaculture model

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

# Stocking conditions: L = 40mm ~ W = 0.4g
  nstart.p = c(P = 5000, L = 40)
  t.p = seq(0, 365*2, 1)
  area = 1000
  