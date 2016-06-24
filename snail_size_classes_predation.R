require(deSolve)

snail_size=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
    S1=n[3]
    S2=n[4]
    S3=n[5]
    
    Bp = (a*(L/10)^b)/10 #Mean prawn mass using conversion from length (cm) to weight (g)
    
    Bn = S1*0.01 + S2*0.1 + S3*0.5
    
      rps1 = Bp / (S1*0.01) #Prawn to snail mass ratio of size class 1
      rps2 = Bp / (S2*0.1)  #Prawn to snail mass ratio of size class 2
      rps3 = Bp / (S3*0.5)  #Prawn to snail mass ratio of size class 3
      
      rps = Bp/Bn #total prawn to snail mass ratio
      
        Th = 0.53623*exp(-0.05167*rps) #from fit to Th / rps data in sokolow 2013
        alpha1 = 1+0.03*rps1
        alpha2 = 1+0.03*rps2
        alpha3 = 1+0.03*rps3
    
    Bm = P*Bp #per prawn biomass conversion to total biomass
    
    dPdt= -P*(muP*L^-0.25 + Bm/phi) #Number of prawns subject to baseline mortality rate and density dependent mortality

    dLdt= k/(1+gam*Bm)*(linf - L) #Mean prawn length growing at growth rate k limited by max length linf, and pop density (B/phi)
    
      #psi1 = (P*alpha*gam1*S1) / (1+alpha*Th*(gam1*S1 + gam2*S2 + gam3*S3))
        psi1 = (P*alpha1) / (1+alpha1*Th*(S1+S2+S3))
      #psi2 = (P*alpha*gam2*S2) / (1+alpha*Th*(gam1*S1 + gam2*S2 + gam3*S3))
        psi2 = (P*alpha2) / (1+alpha2*Th*(S1+S2+S3))
      #psi3 = (P*alpha*gam3*S3) / (1+alpha*Th*(gam1*S1 + gam2*S2 + gam3*S3))
        psi3 = (P*alpha3) / (1+alpha3*Th*(S1+S2+S3))
      
      N = S1+S2+S3
          
    dS1dt = f*(1-N/Kn)*(S2+S3) - muN*S1 - psi1*S1 - g1*S1
    
    dS2dt = g1*S1 - muN*S2 - psi2*S2 - g2*S2
    
    dS3dt = g2*S2 - muN*S3 - psi3*S3
      
    
    return(list(c(P, L, S1, S2, S3)))
  }) 
} 

#Set initial values and parameters ##################
nstart=c(P = 15000, L = 25, S1 = 1000, S2 = 400, S3 = 100)
time=seq(0,365*2,1)

#List parameters and values
parameters=c(
  #PRawn model parameters
    a = 0.096868,
    b = 3.2944,
    gam = 1e-6, #crowding parameter that reduces growth rate at high density
    muP = 0.006136986, #baseline prawn mortality rate
    phi = 40000000, #biomass-assessed density dependence parameter in one hactare
    k = 0.00339726, #growth rate (mm/day)
    linf = 206, #max length (mm)
    
  #Snail model parameters
    f = 0.16, #Birth rate of adult snails
    muN = 1/365, #Natural mortality rate of snails (assume mean lifespan =100 days)
    g1 = 1/37, #Time to transition to 8mm class assuming growth rate of 1.5mm/14 days (McCreesh 2014) 
    g2 = 1/62, #Time to transition to 12mm class assuming growth rate of 0.9mm/14 days (McCreesh 2014) 
    Kn = 15000
  )


#Run & plot ############
output_ss=as.data.frame(ode(nstart,time,snail_size,parameters))

#plot results ###############
  plot(output_ss$time, output_ss$P, type = 'l', col = 'red', lwd=2)
  plot(output_ss$time, output_ss$L, type = 'l', col = 'green', lwd=2)