## Epidemiological model including snail size classes and immigration of snails; no predation

require(deSolve)

#mating function
fx<-function(x, mean.worm, clump){
  alpha<-(mean.worm)/(clump+mean.worm)
  (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
}

phi_Wk<-function(W, phi){
  alpha<-W/(W+phi)
  1-( (1-alpha)^(phi+1) * (integrate(fx, 0, 2*pi, W, phi, stop.on.error = F)$value) /(2*pi)  )
}

#Model with no diagonal transitions, all horizontal and vertical transitions #####
snail_epi_allvh_imm = function(t, n, parameters) { 
  with(as.list(parameters),{
    
    S1=n[1]
    S2=n[2]
    S3=n[3]
    E1=n[4]
    E2=n[5]
    E3=n[6]
    I1=n[7]
    I2=n[8]
    I3=n[9]
    Wt=n[10]
    Wu=n[11]
    
    W = cov*Wt + (1-cov)*Wu
    
    mate = phi_Wk(W = W, phi = phi)  #Mating probability
    S = S1+S2+S3    # Total susceptible snails
    E = E1+E2+E3    # Total exposed snails
    I = I1+I2+I3       # Total infected snails
    N = S+E+I       # Total number of snails
    M = m*0.5*W*mate # Miracidial density per person as a function of mean worm burden (W), miracidial shedding rate (m), and mating probability
    # rho = g2/sigma  # Fraction of E2 snails that transition directly to I3
    
    dS1dt = iota*siteS1 + f*(1-N/Kn)*(S2+S3+z*(E2+E3)) - muN1*S1 - psi1*S1 - g1*S1 - (beta)*M*H*S1 - iota*S1
    
    dS2dt = iota*siteS2 + g1*S1 - muN2*S2 - psi2*S2 - g2*S2 - (beta)*M*H*S2 - iota*S2
    
    dS3dt = iota*siteS3 + g2*S2 - muN3*S3 - psi3*S3 - (beta)*M*H*S3 - iota*S3
    
    dE1dt = iota*siteE1 + (beta)*M*H*S1 - muN1*E1 - psi1*E1 - sigma*E1 - g1*E1 - iota*E1
    
    dE2dt = iota*siteE2 + (beta)*M*H*S2 + g1*E1 - muN2*E2 - psi2*E2 - sigma*E2 - g2*E2 - iota*E2
    
    dE3dt = iota*siteE3 + (beta)*M*H*S3 + g2*E2 - muN3*E3 - psi3*E3 - sigma*E3 - iota*E3
    
    dI1dt = iota*siteI1 + sigma*E1 - (muN1+muI)*I1 - psi2*I1 - g1*I1 - iota*I1

    dI2dt = iota*siteI2 + sigma*E2 + g1*I1 - (muN2+muI)*I2 - psi2*I2 - g2*I2 - iota*I2
    
    dI3dt = iota*siteI3 + sigma*E3 + g2*I2 - (muN3+muI)*I3 - psi3*I3 - iota*I3
    
    dWtdt = lambda*I1 + theta1*lambda*I2 + theta2*lambda*I3 - (muW + muH)*Wt
    
    dWudt = lambda*I1 + theta1*lambda*I2 + theta2*lambda*I3 - (muW + muH)*Wu
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI1dt, dI2dt, dI3dt, dWtdt, dWudt)))
  }) 
} 
