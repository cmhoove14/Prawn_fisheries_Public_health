predation<-function(P, P.L, N1, N2, N3){
  
  N = N1+N2+N3

  P.bm = (0.096868*(P.L/10)^3.2944)/10 #prawn length to weight conversion from Lalrisanga et al (rosenbergii)
  
  N1.bm = 0.3124*0.4^2 - 0.1205*0.4 #4mm size class; length to weight from Sanna's data
  N2.bm = 0.3124*0.8^2 - 0.1205*0.8 #8mm size class; length to weight from Sanna's data
  N3.bm = 0.3124*1.2^2 - 0.1205*1.2 #12mm size class; length to weight from Sanna's data
  
  N.bm = N1.bm*N1 + N2.bm*N2 + N3.bm*N3
    
  rps1 = P.bm / N1.bm
  rps2 = P.bm / N2.bm
  rps3 = P.bm / N3.bm
  
  alpha1 = 1+0.03*rps1 #attack rate as a function of prawn to snail biomass ratio
  alpha2 = 1+0.03*rps2 #attack rate as a function of prawn to snail biomass ratio
  alpha3 = 1+0.03*rps3 #attack rate as a function of prawn to snail biomass ratio
  
  handle1 = (0.53623*exp(-0.05167*rps1))^-1 #handling time as a function of prawn to snail mass ratio
  handle2 = (0.53623*exp(-0.05167*rps2))^-1 #handling time as a function of prawn to snail mass ratio
  handle3 = (0.53623*exp(-0.05167*rps3))^-1 #handling time as a function of prawn to snail mass ratio
  
  gam1 = P.bm / 100
  gam2 = 1/3
  gam3 = (1-1/3)/P.bm + 1/3
    
  gam_star1 = gam1 / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  gam_star2 = gam2 / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  gam_star3 = gam3 / (gam1*(N1/N) + gam2*(N2/N) + gam3*(N3/N))
  
  psi1 = (P*alpha1*gam_star1*N1) / (1+sum((alpha1/handle1)*(gam_star1*N1*(handle1/(handle1+handle2+handle3))),
                                          (alpha2/handle2)*(gam_star2*N2*(handle2/(handle1+handle2+handle3))),
                                          (alpha3/handle3)*(gam_star3*N3)))
  
  psi2 = (P*alpha2*gam_star2*N2) / (1+sum((alpha1/handle1)*(gam_star1*N1*(handle1/(handle1+handle2+handle3))),
                                          (alpha2/handle2)*(gam_star2*N2*(handle2/(handle1+handle2+handle3))),
                                          (alpha3/handle3)*(gam_star3*N3)))
  
  psi3 = (P*alpha3*gam_star3*N3) / (1+sum((alpha1/handle1)*(gam_star1*N1*(handle1/(handle1+handle2+handle3))),
                                          (alpha2/handle2)*(gam_star2*N2*(handle2/(handle1+handle2+handle3))),
                                          (alpha3/handle3)*(gam_star3*N3)))
  
  return(c(psi1,psi))
}