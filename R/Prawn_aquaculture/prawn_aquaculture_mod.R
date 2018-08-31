## Modeling prawn growth and mortality in aquaculture

require(deSolve)

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