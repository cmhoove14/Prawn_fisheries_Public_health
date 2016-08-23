#Sensitivity of snail epi model (reduced to exclude prawn dynamics)

require(sensitivity)
require(ggplot2)

source('snail_epi.R')

#Get parameter sets to sample from #######################
sims3 = 100
#snail parameters
  f.range<-seq(0.05, 0.45, length.out = sims3)
  Kn.range<-seq(20, 60, length.out = sims3)
  z.range<-seq(0, 1, length.out = sims3)
  g1.range<-seq(parameters['g1']/1.2, parameters['g1']*1.2, length.out = sims3)
  g2.range<-seq(parameters['g2']/1.2, parameters['g2']*1.2, length.out = sims3)
  muN1.range<-seq(parameters['muN1']/1.2, parameters['muN1']*1.2, length.out = sims3)
  muN2.range<-seq(parameters['muN2']/1.2, parameters['muN2']*1.2, length.out = sims3)
  muN3.range<-seq(parameters['muN3']/1.2, parameters['muN3']*1.2, length.out = sims3)
  muI.range<-seq(parameters['muI']/1.2, parameters['muI']*1.2, length.out = sims3)
  beta.range<-seq(parameters['beta']/10, parameters['beta']*10, length.out = sims3)
  lambda.range<-seq(parameters['lambda']/10, parameters['lambda']*10, length.out = sims3)
  sigma.range<-seq(1/25, 1/75, length.out = sims3)
  theta.range<-seq(2, 10, length.out = sims3)
  m.range<-seq(parameters['m']/10, parameters['m']*10, length.out = sims3)
  
  paranges3<-cbind(f=f.range, 
                   Kn=Kn.range, 
                   z=z.range, 
                   g1=g1.range, 
                   g2=g2.range, 
                   muN1=muN1.range,
                   muN2=muN2.range,
                   muN3=muN3.range,
                   muI=muI.range,
                   beta=beta.range,
                   lambda=lambda.range,
                   sigma=sigma.range,
                   theta=theta.range,
                   m=m.range)

constantparams<-matrix(ncol = length(parameters), nrow = sims3)

for(i in 1:length(parameters)){
  constantparams[,i] = rep(parameters[i],sims3)
}

colnames(constantparams)<-names(parameters)
vars3<-colnames(paranges3)

#First check scatter plots of outcomes across each parameter range ############
  #i.e. check for monotinicity
  outputfill1<-matrix(ncol = length(vars3), nrow = sims3)
  outputfill2<-matrix(ncol = length(vars3), nrow = sims3)
  outputfill3<-matrix(ncol = length(vars3), nrow = sims3)
  
  nstart3 = c(S1 = snailepiq$S1, S2 = snailepiq$S2, S3 = snailepiq$S3, 
              E1 = snailepiq$E1, E2 = snailepiq$E2, E3 = snailepiq$E3, 
              I2 = snailepiq$I2, I3 = snailepiq$I3, W = snailepiq$W)
  
  for(j in 1:length(vars3)){
    for(i in 1:sims3){
      print(c(j, i))

      other<-constantparams[, -which(colnames(constantparams) %in% vars3[j])]
      test<-paranges3[, which(colnames(paranges3) %in% vars3[j])]
      
      parametersuse<-cbind(other, test) 
      colnames(parametersuse)[dim(parametersuse)[2]]<-vars3[j]
      
      parameters3<-parametersuse[i,]
      
      outputest = as.data.frame(ode(nstart3, time, snail_epi, parameters3))
        
        #par(opar)
        #plot(outputest$time, (outputest$S1 + outputest$S2 + outputest$S3 +
         #                     outputest$E1 + outputest$E2 + outputest$E3 +
          #                    outputest$I2 + outputest$I3), type = 'l', lwd=2,
           #  ylab = 'total snail pop', xlab = 'time (days)', ylim = c(0,60))
      
      outputfill1[i,j] = sum(outputest$S1[dim(outputest)[1]] + outputest$S2[dim(outputest)[1]] + outputest$S3[dim(outputest)[1]] + 
                           outputest$E1[dim(outputest)[1]] + outputest$E2[dim(outputest)[1]] + outputest$E3[dim(outputest)[1]] + 
                           outputest$I2[dim(outputest)[1]] + outputest$I3[dim(outputest)[1]])
      outputfill2[i,j] = sum(outputest$I2[dim(outputest)[1]] + outputest$I3[dim(outputest)[1]])
      outputfill3[i,j] = outputest$W[dim(outputest)[1]]
      
      
    }
    par(mfrow = c(3,1), mar = c(4,3.75,1,0.4)+0.1)
      
      plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill1[,j], 
           xlab = vars3[j], ylab = 'snail pop @10 yrs',
           pch = 16, cex = 0.75, col = 'blue', 
           ylim = c(0,60))
      
      plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill2[,j], 
           xlab = vars3[j], ylab = 'infecteds @10 yrs',
           pch = 16, cex = 0.75, col = 'red', 
           ylim = c(0,3))
      
      plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill3[,j], 
           xlab = vars3[j], ylab = 'mean worm burden @10 yrs',
           pch = 16, cex = 0.75, col = 'purple', 
           ylim = c(0,120))
  }
  
#Update parameter ranges ########
  #reprodction rate was too low
    f.range<-seq(0.1, max(f.range), length.out = sims3)
  #beta was all over the place but appears as though range below might be most feasible
    beta.range<-seq(6e-7, 5e-6, length.out = sims3)
  #m was also all over the place, try the range below
    m.range<-seq(0.25, 2, length.out = sims3)  
  #lambda was close to monotonic but had a slight rise at the beginning, so change its lower value
    lambda.range<-seq(0.012, 0.2, length.out = sims3)
    
#Re-run code above to update monotinicity plots and re-verify monotinicity before running code below ##############
    paranges3<-cbind(f=f.range, 
                     Kn=Kn.range, 
                     z=z.range, 
                     g1=g1.range, 
                     g2=g2.range, 
                     muN1=muN1.range,
                     muN2=muN2.range,
                     muN3=muN3.range,
                     muI=muI.range,
                     beta=beta.range,
                     lambda=lambda.range,
                     sigma=sigma.range,
                     theta=theta.range,
                     m=m.range)
    
    outputfill1<-matrix(ncol = length(vars3), nrow = sims3)
    outputfill2<-matrix(ncol = length(vars3), nrow = sims3)
    outputfill3<-matrix(ncol = length(vars3), nrow = sims3)
    
    nstart3 = c(S1 = snailepiq$S1, S2 = snailepiq$S2, S3 = snailepiq$S3, 
                E1 = snailepiq$E1, E2 = snailepiq$E2, E3 = snailepiq$E3, 
                I2 = snailepiq$I2, I3 = snailepiq$I3, W = snailepiq$W)
#just run parameters that were updated    
    for(j in c(1,10,11,14)){
      for(i in 1:sims3){
        print(c(j, i))
        
        other<-constantparams[, -which(colnames(constantparams) %in% vars3[j])]
        test<-paranges3[, which(colnames(paranges3) %in% vars3[j])]
        
        parametersuse<-cbind(other, test) 
        colnames(parametersuse)[dim(parametersuse)[2]]<-vars3[j]
        
        parameters3<-parametersuse[i,]
        
        outputest = as.data.frame(ode(nstart3, time, snail_epi, parameters3))
        
        #par(opar)
        #plot(outputest$time, (outputest$S1 + outputest$S2 + outputest$S3 +
        #                     outputest$E1 + outputest$E2 + outputest$E3 +
        #                    outputest$I2 + outputest$I3), type = 'l', lwd=2,
        #  ylab = 'total snail pop', xlab = 'time (days)', ylim = c(0,60))
        
        outputfill1[i,j] = sum(outputest$S1[dim(outputest)[1]] + outputest$S2[dim(outputest)[1]] + outputest$S3[dim(outputest)[1]] + 
                                 outputest$E1[dim(outputest)[1]] + outputest$E2[dim(outputest)[1]] + outputest$E3[dim(outputest)[1]] + 
                                 outputest$I2[dim(outputest)[1]] + outputest$I3[dim(outputest)[1]])
        outputfill2[i,j] = sum(outputest$I2[dim(outputest)[1]] + outputest$I3[dim(outputest)[1]])
        outputfill3[i,j] = outputest$W[dim(outputest)[1]]
        
        
      }
      par(mfrow = c(3,1), mar = c(4,3.75,1,0.4)+0.1)
      
      plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill1[,j], 
           xlab = vars3[j], ylab = 'snail pop @10 yrs',
           pch = 16, cex = 0.75, col = 'blue', 
           ylim = c(0,60))
      
      plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill2[,j], 
           xlab = vars3[j], ylab = 'infecteds @10 yrs',
           pch = 16, cex = 0.75, col = 'red', 
           ylim = c(0,3))
      
      plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill3[,j], 
           xlab = vars3[j], ylab = 'mean worm burden @10 yrs',
           pch = 16, cex = 0.75, col = 'purple', 
           ylim = c(0,120))
    }
    
#Augment LHc ##########
#create a matrix of indices for the LHC where nrow=number of sims and ncol=number of variables
  LHC_indices3<-matrix(0, nrow=sims3, ncol=14)  
#Add structure of latin hypercube
  for(j in 1:dim(LHC_indices3)[2]){
    LHC_indices3[,j]<-sample(1:sims3, size=sims3, replace=FALSE)
  }  
#Fill LHC values
  params03<-cbind(f=f.range[LHC_indices3[,1]], 
                  Kn=Kn.range[LHC_indices3[,2]], 
                  z=z.range[LHC_indices3[,3]], 
                  g1=g1.range[LHC_indices3[,4]], 
                  g2=g2.range[LHC_indices3[,5]], 
                  muN1=muN1.range[LHC_indices3[,6]],
                  muN2=muN2.range[LHC_indices3[,7]],
                  muN3=muN3.range[LHC_indices3[,8]],
                  muI=muI.range[LHC_indices3[,9]],
                  beta=beta.range[LHC_indices3[,10]],
                  lambda=lambda.range[LHC_indices3[,11]],
                  sigma=sigma.range[LHC_indices3[,12]],
                  theta=theta.range[LHC_indices3[,13]],
                  m=m.range[LHC_indices3[,14]])
  
  params3<-constantparams[, -which(colnames(constantparams) %in% vars3)]
  
  #Add in sampled values of parameters of interest from above
  params.fin3<-cbind(params3, params03)  
  
#Run model through with parameter sets, save outcome variables and merge with params ##################
  snails<-numeric()
  infs<-numeric()
  ws<-numeric()

  time3<-seq(0, 365*20, by=1)
  nstart3= c(S1 = 10*area, 
             S2 = 0, 
             S3 = 0, 
             E1 = 0, 
             E2 = 0, 
             E3 = 0, 
             I2 = 0, 
             I3 = 0, 
             W = 2)
  
  for(i in 1:sims3){
    print(i)
    parameters3<-params.fin3[i,]
    
    output3 = as.data.frame(ode(nstart3, time3, snail_epi, parameters3))
    
    plot(output3$time, (output3$S1 + output3$S2 + output3$S3 + 
                        output3$E1 + output3$E2 + output3$E3 + 
                        output3$I2 + output3$I3), type = 'l', lwd=2, xlab = 'time', ylab = 'snails')
      lines(output3$time, (output3$I2 + output3$I3), col = 'red', lwd = 2)
    
    snails[i] = sum(output3$S1[dim(output3)[1]], output3$S2[dim(output3)[1]], output3$S3[dim(output3)[1]], 
                    output3$E1[dim(output3)[1]], output3$E2[dim(output3)[1]], output3$E3[dim(output3)[1]], 
                    output3$I2[dim(output3)[1]], output3$I3[dim(output3)[1]])
    infs[i] = sum(output3$I2[dim(output3)[1]], output3$I3[dim(output3)[1]])
    ws[i] = output3$W[dim(output3)[1]] 
    print(c(snails[i], infs[i], ws[i]))
  }  
  
  
  params.snail<-as.data.frame(cbind(params.fin3, snails, infs, ws))
  snail.test<-cbind(params.snail[, -which(colnames(params.snail) %in% colnames(params3))])
  
  #save default plot settings
  opar<-par()
  
#Check scatter plots to assess sensitivity ############## 
  #of harvest size and time, mean length at harvest, and pop size to tested parameters    
  par(mfrow = c(2,7), mar = c(4,3.75,1,0.4)+0.1)
  
  for(i in 1:length(vars3)){
    plot(x = snail.test[,i], y = snail.test[,dim(snail.test)[2]], pch = 16, col = 'black', cex = 0.7,
         xlab = vars3[i], ylab = 'Snail population')
  }
  
  for(i in 1:length(vars3)){
    plot(x = snail.test[,i], y = snail.test[,dim(snail.test)[2]-1], pch = 16, col = 'red', cex = 0.7,
         xlab = vars3[i], ylab = 'Infected snails')
  }
  
  for(i in 1:length(vars3)){
    plot(x = snail.test[,i], y = snail.test[,dim(snail.test)[2]-2], pch = 16, col = 'purple', cex = 0.7,
         xlab = vars3[i], ylab = 'Mean worm burden')
  }
