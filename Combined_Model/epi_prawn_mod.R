source('Epi_Model/snail_epi_mod.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')

snail_prawn_model = function(t, n, parameters) {  
  with(as.list(parameters),{
    
    S1=n[1]
    S2=n[2]
    S3=n[3]
    E1=n[4]
    E2=n[5]
    E3=n[6]
    I2=n[7]
    I3=n[8]
    W=n[9]
    P=n[10]
    L=n[11]
    
    N1 = S1+E1     # Total snails of size class 1
    N2 = S2+E2+I2  # Total snails of size class 2
    N3 = S3+E3+I3  # Total snails of size class 3
    N = N1+N2+N3   # Total number of snails
    
    # Miracidial density per person as a function of mean worm burden (W) and miracidial shedding rate (m)
    mate = phi_Wk(W = W, k = k)  #Mating probability
    M = m*0.5*W*mate
    
    # Mean and total prawn biomass, converting from length (mm) to weight (g)
    Bm.p = (a.p/10)*(L/10)^b.p
    Bm.t = P*Bm.p
    
    # Mean snail biomass in each size class, converting from length (cm) to weight (g)
    Bm.n1 = a.s*0.4^b.s  # 4mm size class
    Bm.n2 = a.s*0.8^b.s  # 8mm size class
    Bm.n3 = a.s*1.2^b.s  # 12mm size class
    
    # Prawn-to-snail mass ratios for each size class
    Bm.r1 = Bm.p / Bm.n1
    Bm.r2 = Bm.p / Bm.n2
    Bm.r3 = Bm.p / Bm.n3
    
    # Attack rates for each size class as a function of biomass ratio (4.5 is biomass ratio below which refuge exists)
    alpha1 = ifelse(Bm.r1 > 3, ar.slp*log(Bm.r1), 0)  #if biomass ratio <3, size refuge exists (from Sokolow 2014)
    alpha2 = ifelse(Bm.r2 > 3, ar.slp*log(Bm.r2), 0)  #if biomass ratio <3, size refuge exists (from Sokolow 2014)
    alpha3 = ifelse(Bm.r3 > 3, ar.slp*log(Bm.r3), 0)  #if biomass ratio <3, size refuge exists (from Sokolow 2014)
    
    #alpha1 = ifelse(-log(3) + ar*log(Bm.r1) < 0, -log(3) + ar*log(Bm.r1), 0)
    
    # Adjusted attack rates, accounting for area of interest and limiting factors in the wild
    alpha_star1 = alpha1/sqrt(A)
    alpha_star2 = alpha2/sqrt(A)
    alpha_star3 = alpha3/sqrt(A)
    
    # Handling times for each size class as a function of biomass ratio
    handle1 = 1/(th*Bm.r1)
    handle2 = 1/(th*Bm.r2)
    handle3 = 1/(th*Bm.r3)
    
    # Functional responses for each size/infection class
    psiS1 = (alpha_star1*(S1/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiS2 = (alpha_star2*(S2/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiS3 = (alpha_star3*(S3/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiE1 = (alpha_star1*(E1/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiE2 = (alpha_star2*(E2/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiE3 = (alpha_star3*(E3/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiI2 = (alpha_star2*(I2/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    psiI3 = (alpha_star3*(I3/A)) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A)))
    
    #psi1 = (P*alpha_star1) / (1+sum(alpha_star1*handle1*N1, alpha_star2*handle2*N2, alpha_star3*handle3*N3))
    #psi2 = (P*alpha_star2) / (1+sum(alpha_star1*handle1*N1, alpha_star2*handle2*N2, alpha_star3*handle3*N3))
    #psi3 = (P*alpha_star3) / (1+sum(alpha_star1*handle1*N1, alpha_star2*handle2*N2, alpha_star3*handle3*N3))
    
    # Fraction of E2 snails that transition directly to I3
    rho = g2/sigma
    
    ## Model equations:
    
    # Susceptible snails, class 1
    dS1dt = f*(1-N/Kn)*(S2+S3+z*(E2+E3)) - muN1*S1 - psiS1*P - g1*S1 - beta*M*H*S1
    
    # Susceptible snails, class 2
    dS2dt = g1*S1 - muN2*S2 - psiS2*P - g2*S2 - beta*M*H*S2
    
    # Susceptible snails, class 3
    dS3dt = g2*S2 - muN3*S3 - psiS3*P - beta*M*H*S3
    
    # Exposed (prepatent) snails, class 1
    dE1dt = beta*M*H*S1 - muN1*E1 - psiE1*P - sigma*E1
    
    # Exposed (prepatent) snails, class 2
    dE2dt = beta*M*H*S2 - muN2*E2 - psiE2*P - sigma*E2
    
    # Exposed (prepatent) snails, class 3
    dE3dt = beta*M*H*S3 - muN3*E3 - psiE3*P - sigma*E3
    
    # Infectious (shedding) snails, class 2
    dI2dt = sigma*E1 + sigma*(1-rho)*E2 - (muN2 + muI)*I2 - psiI2*P
    
    # Infectious (shedding) snails, class 3
    dI3dt = sigma*rho*E2 + sigma*E3 - (muN3 + muI)*I3 - psiI3*P
    
    # Mean human worm burden
    dWdt = lambda*I2/A + theta*lambda*I3/A - (muW + muH)*W
    
    # Prawn abundance
    dPdt = -P*(muP*Bm.p^d + phi*Bm.t)
    
    # Mean prawn length (mm)
    dLdt = k/(1+gam*Bm.t)*(linf - L)
    
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI2dt, dI3dt, dWdt, dPdt, dLdt)))
  }) 
} 
