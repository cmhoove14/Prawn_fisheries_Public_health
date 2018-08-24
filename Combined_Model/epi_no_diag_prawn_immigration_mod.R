source('Epi_Model/snail_epi_mod_no_diag_immigration.R')
source('Prawn_aquaculture/prawn_aquaculture_mod.R')

snail_prawn_model_imm = function(t, n, parameters) {  
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
    P=n[12]
    L=n[13]
    
    N1 = S1+E1     # Total snails of size class 1
    N2 = S2+E2+I2  # Total snails of size class 2
    N3 = S3+E3+I3  # Total snails of size class 3
    N = N1+N2+N3   # Total number of snails
    W = cov*Wt + (1-cov)*Wu  #Mean worm burden weighted between treated and untreated populations
    
    # Miracidial density per person as a function of mean worm burden (W) and miracidial shedding rate (m)
    mate = phi_Wk(W = W, phi = phi)  #Mating probability
    M = m*0.5*W*mate
    
    # Mean and total prawn biomass, converting from length (cm) to weight (g)
    Bm.p = 10^a.p*(L/10)^b.p
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
    alpha_star1 = alpha1/eps
    alpha_star2 = alpha2/eps
    alpha_star3 = alpha3/eps
    
    # Handling times for each size class as a function of biomass ratio
    handle1 = 1/(th*Bm.r1)
    handle2 = 1/(th*Bm.r2)
    handle3 = 1/(th*Bm.r3)
    
    # Functional responses for each size/infection class
    psiS1 = (alpha_star1*(S1/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiS2 = (alpha_star2*(S2/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiS3 = (alpha_star3*(S3/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiE1 = (alpha_star1*(E1/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiE2 = (alpha_star2*(E2/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiE3 = (alpha_star3*(E3/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiI1 = (alpha_star1*(I1/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiI2 = (alpha_star2*(I2/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    psiI3 = (alpha_star3*(I3/A)^nfr) / (1 + sum(alpha_star1*handle1*(S1/A),
                                            alpha_star2*handle2*(S2/A),
                                            alpha_star3*handle3*(S3/A),
                                            alpha_star1*handle1*(E1/A),
                                            alpha_star2*handle2*(E2/A),
                                            alpha_star3*handle3*(E3/A),
                                            alpha_star1*handle1*(I1/A),
                                            alpha_star2*handle2*(I2/A),
                                            alpha_star3*handle3*(I3/A))^nfr)
    
    #print(c(psiS1, psiS2, psiS3, psiE1, psiE2, psiE3, psiI1, psiI2, psiI3)*P)
    
    ## Model equations:
    dS1dt = iota*siteS1 + f*(1-N/Kn)*(S2+S3+z*(E2+E3)) - muN1*S1 - psiS1*P - g1*S1 - (beta)*M*H*S1 - iota*S1
    
    dS2dt = iota*siteS2 + g1*S1 - muN2*S2 - psiS2*P - g2*S2 - (beta)*M*H*S2 - iota*S2
    
    dS3dt = iota*siteS3 + g2*S2 - muN3*S3 - psiS3*P - (beta)*M*H*S3 - iota*S3
    
    dE1dt = iota*siteE1 + (beta)*M*H*S1 - muN1*E1 - psiE1*P - sigma*E1 - g1*E1 - iota*E1
    
    dE2dt = iota*siteE2 + (beta)*M*H*S2 + g1*E1 - muN2*E2 - psiE2*P - sigma*E2 - g2*E2 - iota*E2
    
    dE3dt = iota*siteE3 + (beta)*M*H*S3 + g2*E2 - muN3*E3 - psiE3*P - sigma*E3 - iota*E3
    
    dI1dt = iota*siteI1 + sigma*E1 - (muN1+muI)*I1 - psiI1*P - g1*I1 - iota*I1

    dI2dt = iota*siteI2 + sigma*E2 + g1*I1 - (muN2+muI)*I2 - psiI2*P - g2*I2 - iota*I2
    
    dI3dt = iota*siteI3 + sigma*E3 + g2*I2 - (muN3+muI)*I3 - psiI3*P - iota*I3
    
    dWtdt = lambda*I1/A + theta1*lambda*I2/A + theta2*lambda*I3/A - (muW + muH)*Wt
    
    dWudt = lambda*I1/A + theta1*lambda*I2/A + theta2*lambda*I3/A - (muW + muH)*Wu

    # Prawn abundance
    dPdt = -P*(muP*Bm.p^d + om*Bm.t)
    
    # Mean prawn length (mm)
    dLdt = k/(1+gam*Bm.t)*(linf - L)
    
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI1dt, dI2dt, dI3dt, dWtdt, dWudt, dPdt, dLdt)))
  }) 
} 
