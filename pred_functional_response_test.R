## Modeling functional response of prawn on different size classes of snails

## To-do list:
  # Continue to tweak attack rate and handling time functions (i.e. to capture size refuge)
  # Resolve scale issue

predation = function(P, P.L, N1, N2, N3) {
  
  N = N1 + N2 + N3
  
  P.bm = (0.096868*(P.L/10)^3.2944)/10 # Prawn length (mm) to weight (g) conversion,
                                       # based on Lalrisanga et al. 2012 (for M. rosenbergii)
  
  N1.bm = 0.187178454*0.4^2.536764792 # 4mm size class; length-to-weight conversion based on fitting allometric equation to Sanna's data on B. glabrata
  N2.bm = 0.187178454*0.8^2.536764792 # 8mm size class
  N3.bm = 0.187178454*1.2^2.536764792 # 12mm size class
                                      # Alternative polynomial function: 0.3124*x^2 - 0.1205*x
  
  # Ratios of prawn-to-snail biomass
  rps1 = P.bm / N1.bm
  rps2 = P.bm / N2.bm
  rps3 = P.bm / N3.bm
  
  # Attack rates as a function of biomass ratio
  alpha1 = ifelse(-log(3) + 1.1906*log(rps1) > 0, -log(3) + 1.1906*log(rps1), 0) # Alternative linear function: 0.037192*rps
  alpha2 = ifelse(-log(3) + 1.1906*log(rps2) > 0, -log(3) + 1.1906*log(rps2), 0)
  alpha3 = ifelse(-log(3) + 1.1906*log(rps3) > 0, -log(3) + 1.1906*log(rps3), 0)
  
  # Handling times as a function of biomass ratio
  handle1 = 1/(0.38561*rps1) # Alternative function: 0.53623*exp(-0.05167*rps)
  handle2 = 1/(0.38561*rps2)
  handle3 = 1/(0.38561*rps3)
  
  # Preference variables, currently ignoring these
  gam1 = 1 #1/3 + 0.5*exp(-0.15*P.bm)
  gam3 = 1 #1/3*exp(-12*exp(-0.2*P.bm))
  gam2 = 1 #1 - gam1 - gam3
  
  # Adjusted preference variables, currently ignoring these
  # gam_star1 = gam1*(N1/N) / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  # gam_star2 = gam2*(N2/N) / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  # gam_star3 = gam3*(N3/N) / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  
  # Size-class specific functional responses
  psi1 = (P*alpha1*gam1*N1) / (1+sum(alpha1*handle1*gam1*N1,
                                     alpha2*handle2*gam2*N2,
                                     alpha3*handle3*gam3*N3))
  
  psi2 = (P*alpha2*gam2*N2) / (1+sum(alpha1*handle1*gam1*N1,
                                     alpha2*handle2*gam2*N2,
                                     alpha3*handle3*gam3*N3))
  
  psi3 = (P*alpha3*gam3*N3) / (1+sum(alpha1*handle1*gam1*N1,
                                     alpha2*handle2*gam2*N2,
                                     alpha3*handle3*gam3*N3))
  
  # Aggregate functional response
  psi = psi1+psi2+psi3
  
  return(list("psi1"=psi1, "psi2"=psi2, "psi3"=psi3, "psi"=psi))
}

# Plot functional response for a given prawn size and snail population
x = c(0:50)
fr = matrix(ncol = 4, nrow = length(x))
for (i in x) {
  pred = predation(P = 1, P.L = 75, N1 = 0.4*i, N2 = 0.3*i, N3 = 0.3*i)
  fr[i+1,1] = pred$psi
  fr[i+1,2] = pred$psi1
  fr[i+1,3] = pred$psi2
  fr[i+1,4] = pred$psi3
  
}
plot(x, fr, type="l", xlab = "Snail density", ylab = "Functional response")
  for(j in 2:4){
    lines(x, fr[,j], col = j)
  }


