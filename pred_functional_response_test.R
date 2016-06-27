predation<-function(P, P.L, N1, N2, N3){
  
  N = N1 + N2 + N3
  
  P.bm = (0.096868*(P.L/10)^3.2944)/10 #prawn length (mm) to weight (g) conversion from Lalrisanga et al (rosenbergii)
  
  N1.bm = 0.3124*0.4^2 - 0.1205*0.4 #4mm size class; length to weight from Sanna's data
  N2.bm = 0.3124*0.8^2 - 0.1205*0.8 #8mm size class; length to weight from Sanna's data
  N3.bm = 0.3124*1.2^2 - 0.1205*1.2 #12mm size class; length to weight from Sanna's data
  
  rps1 = P.bm / N1.bm
  rps2 = P.bm / N2.bm
  rps3 = P.bm / N3.bm
  
  ## TODO: refit functions for alpha and Th
  alpha1 = 1+0.03*rps1 #attack rate as a function of prawn to snail biomass ratio
  alpha2 = 1+0.03*rps2 #attack rate as a function of prawn to snail biomass ratio
  alpha3 = 1+0.03*rps3 #attack rate as a function of prawn to snail biomass ratio
  
  handle1 = (0.53623*exp(-0.05167*rps1))^-1 #handling time as a function of prawn to snail mass ratio
  handle2 = (0.53623*exp(-0.05167*rps2))^-1 #handling time as a function of prawn to snail mass ratio
  handle3 = (0.53623*exp(-0.05167*rps3))^-1 #handling time as a function of prawn to snail mass ratio
  
  gam1 = 1 #1/3 + 0.5*exp(-0.15*P.bm)
  gam3 = 1 #1/3*exp(-12*exp(-0.2*P.bm))
  gam2 = 1 #1 - gam1 - gam3
  
  ## Probably needs tweaking
  #gam_star1 = gam1*(N1/N) / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  #gam_star2 = gam2*(N2/N) / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  #gam_star3 = gam3*(N3/N) / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  
  psi1 = (P*alpha1*gam1*N1) / (1+sum((alpha1/handle1)*gam1*N1,
                                     (alpha2/handle2)*gam2*N2,
                                     (alpha3/handle3)*gam3*N3))
  
  psi2 = (P*alpha2*gam2*N2) / (1+sum((alpha1/handle1)*gam1*N1,
                                     (alpha2/handle2)*gam2*N2,
                                     (alpha3/handle3)*gam3*N3))
  
  psi3 = (P*alpha3*gam3*N3) / (1+sum((alpha1/handle1)*gam1*N1,
                                     (alpha2/handle2)*gam2*N2,
                                     (alpha3/handle3)*gam3*N3))
  
  psi = psi1+psi2+psi3
  
  return(psi)
}

## TODO: experiment with curves for different prawn sizes and snail populations (single or mixed size)
fr = rep.int(0, 51)
for (i in 0:150) {
  fr[i+1] = predation(P = 1, P.L = 75, N1 = i, N2 = 0, N3 = 0)
}
plot(c(0:150), fr, type="l")



