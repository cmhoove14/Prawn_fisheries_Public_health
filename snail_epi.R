#Epidemiological model of different size classes

#To-do list
  #Draw parameter estimates from literature; current ones are taken from Sokolow PNAS paper or from thin air
  #Make sure snail dynamics between each size/infection class are appropriate; 
      #should Infected snails be able to grow?
      #should exposed snails be able to reproduce at a reduced rate?
      #should smallest size class be able to reach infected stage?
  #Come up with R0 expression or some other way to estimate transmission risk

require(deSolve)

snail_epi=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    #Insert prawn density and length equations here once integrating full model
    S1=n[1]
    S2=n[2]
    S3=n[3]
    E1=n[4]
    E2=n[5]
    E3=n[6]
    I2=n[7]
    I3=n[8]
    C=n[9]
    W=n[10]
    
    N = S1+S2+S3+E1+E2+E3+I2+I3 #Total number of snails
    S = S1+S2+S3 #Total susceptible snails
    E = E1+E2+E3 #Total exposed snails
    I = I2+I3 #Total infected snails
    M = m*0.5*W #Miracidial density as a function of mean worm burden (W) and miracidial shedding rate (m); 
                #ADD MATING FUNCTION????
    
    #insert functional responses of predation to each snail size class here
      #psi1 = 
      #psi2 = 
      #psi3 = 
    
    dS1dt = f*(1-N/Kn)*(S2+S3) - muN*S1 - psi1*S1 - g1*S1 - beta*M*S1
    
    dS2dt = g1*S1 - muN*S2 - psi2*S2 - g2*S2 - beta*M*S2
    
    dS3dt = g2*S2 - muN*S3 - psi3*S3 - beta*M*S3
    
    dE1dt = beta*M*S1 - muN*E1 - psi1*E1 - g1*E1
    
    dE2dt = beta*M*S2 + g1*E1 - muN*E2 - psi2*E2 - g2*E2 - sigma*E2
    
    dE3dt = beta*M*S3 + g2*E2 - muN*E3 - psi3*E3 - sigma*E3
    
    dI2dt = sigma*E2 - (muN + muI)*I2 - psi2*I2 - g2*I2
    
    dI3dt = sigma*E3 + g2*I2 - (muN + muI)*I3 - psi3*I3
    
    dCdt = theta2*I2 + theta3*I3 - muC*C
    
    dWdt = lamda*C - (muW+muH)*W
    
    
    return(list(c(dS1dt, dS2dt, dS3dt, dE1dt, dE2dt, dE3dt, dI2dt, dI3dt, dCdt, dWdt)))
  }) 
} 

#Set initial values and parameters ##################
nstart=c(S1 = 1, S2 = 0, S3 = 0, E1 = 0, E2 = 0, E3 = 0, I2 = 0, I3 = 0, C = 0, W = 2)
time=seq(0,365*50,1)

#List parameters and values
parameters=c(
  #Model parameters (in order of appearance)
    m = 1000, #miracidial shedding rate per adult female worm
    f = 0.16, #Birth rate of adult snails (snails/reproductive snail/day; including survival to detection; more like a recruitment rate)
    Kn = 50, #carrying capacity of snails (snails/m^2)
    muN = 1/50, #Natural mortality rate of snails (deaths/snail/day; assume mean lifespan = 100 days)
    psi1 = 0, #TBD: predation rate of prawn cohort on smallest size class of snails
    g1 = 1/37, #Growth rate of smallest size class snails (translated from mm/day from McCreesh 2014 to days of growth required for transition to next size class) 
    beta = 1e-6, #man-to-snail infection probability; (exposed snails/miracidia/snail/day)
    psi2 = 0, #TBD: predation rate of prawn cohort on medium size class of snails
    g2 = 1/62, ##Growth rate of medium size class snails (translated from mm/day from McCreesh 2014 to days of growth required for transition to next size class) 
    psi3 = 0, #TBD: predation rate of prawn cohort on large size class of snails
    sigma = 1/40, #latent period of exposed snails (infected snails/exposed snail/day; assumed to be average of 40 days)
    muI = 1/15, #Additional daily mortality rate of infected snails as a result of infection
    theta2 = 500, #shedding rate of medium infected snails (cercariae/infected medium snail/day)
    theta3 = 1000, #shedding rate of large infected snails (cercariae/infected large snail/day)
    muC = 0.99, #mortality rate of cercariae (dead cercariae/cercariae/day; very small chance of cercariae lasting >1 day)
    lamda = 1e-4, #probability a cercariae interacts with a human and successfully infects them and matures to reproductive age
    muW = 1/(3.3*365), #natural mortality rate of adult worms in humans; assume average lifespan of 3.3 years
    muH = 1/(60*365) #natural mortality of humans that contributes to worm mortality; assume average lifespan of 60 years
)


#Run & plot ##############
output_s.epi=as.data.frame(ode(nstart,time,snail_epi,parameters))
  output_s.epi$S.t = output_s.epi$S1 + output_s.epi$S2 + output_s.epi$S3 #add total susceptible snails to data frame
  output_s.epi$E.t = output_s.epi$E1 + output_s.epi$E2 + output_s.epi$E3 #add total exposed snails to data frame
  output_s.epi$I.t = output_s.epi$I2 + output_s.epi$I3 #add total susceptible snails to data frame
  output_s.epi$N.t = output_s.epi$S.t + output_s.epi$E.t + output_s.epi$I.t #Add total snails between all size & infection classes to data frame
  output_s.epi$t.1 = output_s.epi$S1 + output_s.epi$E1 #add total snails of size class 1 to data frame
  output_s.epi$t.2 = output_s.epi$S2 + output_s.epi$E2 + + output_s.epi$I2 #add total snails of size class 2 to data frame
  output_s.epi$t.3 = output_s.epi$S3+ output_s.epi$E3 + + output_s.epi$I3 #add total snails of size class 2 to data frame
  
par(mfrow = c(2,2))  
plot(x = output_s.epi$time, y = output_s.epi$N.t, type = 'l', col = 'black', lwd=2, xlab = 'time', 
     ylab = 'snail size classes', ylim = c(0,max(output_s.epi$N.t)),
     main = 'Snail Size Classes')
  lines(output_s.epi$time, output_s.epi$t.1, col = 'green', lwd = 2)
  lines(output_s.epi$time, output_s.epi$t.2, col = 'blue', lwd = 2)
  lines(output_s.epi$time, output_s.epi$t.3, col = 'red', lwd = 2)
  #legend('topright', legend = c('total', '1', '2', '3'), lwd = 2, col = c('black', 'green', 'blue', 'red'), cex = 0.7)
  
plot(x = output_s.epi$time, y = output_s.epi$N.t, type = 'l', col = 'black', lwd=2, xlab = 'time', 
       ylab = 'snail infection classes', ylim = c(0,max(output_s.epi$N.t)), 
     main = 'Snail infection classes')
  lines(output_s.epi$time, output_s.epi$I.t, col = 'red', lwd = 2)
  lines(output_s.epi$time, output_s.epi$E.t, col = 'orange', lwd = 2)
  lines(output_s.epi$time, output_s.epi$S.t, col = 'green', lwd = 2)
  #legend('topright', legend = c('total', 'S', 'E', 'I'), lwd = 2, col = c('black', 'green', 'orange', 'red'), cex = 0.7)
  
plot(x = output_s.epi$time, y = output_s.epi$C, type = 'l', col = 'orange', lwd=2, xlab = 'time', 
     ylab = 'cercariae', ylim = c(0,max (output_s.epi$C)),
     main = 'Cercarial density')

plot(x = output_s.epi$time, y = output_s.epi$W, type = 'l', col = 'red', lwd=2, xlab = 'time', 
          ylab = 'worm burden', ylim = c(0,max (output_s.epi$W)),
     main = 'Worm burden')