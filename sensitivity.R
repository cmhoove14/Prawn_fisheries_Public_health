source('snail_prawn_full.R')
require(sensitivity)
require(ggplot2)

#Generate ranges for input parameters of interest ##################
  #Note: these parameter ranges are simply +/- 2 orders of magnitude from the point estimate;
  #could be updated with more informed max/min values and plausible distributions
sims=1000

#snail parameters
  beta.range<-seq(parameters['beta']/100, parameters['beta']*100, length.out = sims)
  lamda.range<-seq(parameters['lambda']/100, parameters['lambda']*100, length.out = sims)
  sigma.range<-seq(parameters['sigma']/100, parameters['sigma']*100, length.out = sims)
  f.range<-seq(parameters['f']/100, parameters['f']*100, length.out = sims)
  z.range<-seq(parameters['z']/100, parameters['z']*100, length.out = sims)
  theta.range<-seq(parameters['theta']/100, parameters['theta']*100, length.out = sims)
#Prawn parameters
  s.range<-seq(parameters['s']/100, parameters['s']*100, length.out = sims)
  phi.range<-seq(parameters['phi']/100, parameters['phi']*100, length.out = sims)
  gam.range<-seq(parameters['gam']/100, parameters['gam']*100, length.out = sims)
  
#create a matrix of indices for the LHC where nrow=number of samples and ncol=number of variables
  LHC_indices<-matrix(0, nrow=sims, ncol=9)  
#Add structure of latin hypercube
  for(j in 1:dim(LHC_indices)[2]){
    LHC_indices[,j]<-sample(1:sims, size=sims, replace=FALSE)
  }  
#Fill LHC values
  params00<-cbind(beta=beta.range[LHC_indices[,1]], 
                  lambda=lamda.range[LHC_indices[,2]], 
                  sigma=sigma.range[LHC_indices[,3]], 
                  f=f.range[LHC_indices[,4]], 
                  z=z.range[LHC_indices[,5]], 
                  theta=theta.range[LHC_indices[,6]], 
                  s=s.range[LHC_indices[,7]], 
                  phi=phi.range[LHC_indices[,8]],
                  gam=gam.range[LHC_indices[,9]])

#Hold all other parameters constant  
  constantparams<-matrix(ncol = length(parameters), nrow = sims)
  
  for(i in 1:length(parameters)){
    constantparams[,i] = rep(parameters[i],sims)
  }
  
  colnames(constantparams)<-names(parameters)
  vars<-colnames(params00)
  params<-constantparams[, -which(colnames(constantparams) %in% vars)]

#Add in sampled values of parameters of interest from above
  params.fin<-cbind(params, params00)
#This provides a data frame to sample from to assess the POIs from above  
  
#First assess sensitivity of harvest time and harvest size ######################
  hs<-numeric()
  ht<-numeric()
  
  area = 10000
  time1 = seq(0, 365*2, 1)
  
  for(i in 1:sims){
    parameters1<-params.fin[i,]
    nstart1= c(S1 = 3.54*area, S2 = 1.02*area, S3 = 0.22*area, 
               E1 = 3.51*area, E2 = 1.1*area, E3 = 0.25*area, 
               I2 = 1.13*area, I3 = 0.23*area, 
               W = 72, 
               P = 5000*area/10000, L = 25)
    
    output = as.data.frame(ode(nstart1, time1, snail_prawn_model, parameters1))
    output$B = ((parameters['a.p']*(output$L/10)^parameters['b.p'])/10)       
    output$Bt = output$B*output$P                                             
    harvest.size = output$B[output$Bt == max(output$Bt)]                      
    harvest.time = output$time[output$Bt == max(output$Bt)]    
    hs[i] = harvest.size
    ht[i] = harvest.time
  }
  
#add harvest time and size to corresponding parameter sets in data frame
  params.harvest<-as.data.frame(cbind(params.fin, hs, ht))
  ph.test<-cbind(params.harvest[, -which(colnames(params.harvest) %in% colnames(params))])
  
  #save default plot settings
    opar<-par()
  
#Check scatter plots to assess sensitivity of harvest size and time to tested parameters    
  par(mfrow = c(3,3), mar = c(4,3.75,1,0.4)+0.1)

  for(i in 1:length(vars)){
      plot(x = ph.test[,i], y = ph.test[,dim(ph.test)[2]], pch = 16, col = 'blue', cex = 0.7,
           xlab = vars[i], ylab = 'Harvest time')
  }
  
  for(i in 1:length(vars)){
    plot(x = ph.test[,i], y = ph.test[,dim(ph.test)[2]-1], pch = 16, col = 'red', cex = 0.7,
          xlab = vars[i], ylab = 'Harvest size')
  }
  
  #restore original plot settings
    par(opar)
    
#Formally analyze sensitivity of harvest time and size using PRCC ############
ht.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
            y = as.data.frame(ph.test[,c(dim(ph.test)[2])]),
            rank = TRUE, conf = 0.95, nboot = 500)

ht.pcc.df<-ht.pcc$PRCC
ht.pcc.df$var = vars
  colnames(ht.pcc.df)[c(1:5)] = c('rank', 'bias', 'std.err', 'min95ci', 'max95ci')

ggplot(ht.pcc.df, aes(x = var, y = rank)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_bar(fill = 'blue', stat = 'identity', width = 0.25)+
  labs(title = 'Harvest time sensitivity', x = 'Parameter', y = 'PRCC') + 
  geom_errorbar(aes(x = var, ymin = min95ci, ymax = max95ci), width = 0.1)
  
    
hs.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
            y = as.data.frame(ph.test[,c(dim(ph.test)[2]-1)]),
            rank = TRUE, conf = 0.95, nboot = 500)

hs.pcc.df<-hs.pcc$PRCC
hs.pcc.df$var = vars
colnames(hs.pcc.df)[c(1:5)] = c('rank', 'bias', 'std.err', 'min95ci', 'max95ci')

ggplot(hs.pcc.df, aes(x = var, y = rank)) +
  theme_bw()+
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
  geom_bar(fill = 'red', stat = 'identity', width = 0.25)+
  labs(title = 'Harvest size sensitivity', x = 'Parameter', y = 'PRCC') + 
  geom_errorbar(aes(x = var, ymin = min95ci, ymax = max95ci), width = 0.1)