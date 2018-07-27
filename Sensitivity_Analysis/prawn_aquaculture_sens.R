source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')
load('Prawn_aquaculture/aquaculture_sims.RData')
require(ggplot2)
require(sensitivity)
require(Rmisc)

#Start with parameters that affect prawn aquaculture dynamics ##################
#  Vary each parameter +/- 25% of it's original value
  sims=1000
  nstart.v = c(P = opt.vol$P_nought, L = 38)
  eta = predict(eta.lm, newdata = data.frame(dens = nstart.v['P']/area))
  
#Prawn parameters
  a.p.range<-seq(par.aqua['a.p']*0.75, par.aqua['a.p']*1.25, length.out = sims)
  b.p.range<-seq(par.aqua['b.p']*0.75, par.aqua['b.p']*1.25, length.out = sims)
  k.range<-seq(par.aqua['k']*0.75, par.aqua['k']*1.25, length.out = sims)
  linf.range<-seq(par.aqua['linf']*0.75, par.aqua['linf']*1.25, length.out = sims)
  muP.range<-seq(par.aqua['muP']*0.75, par.aqua['muP']*1.25, length.out = sims)
  d.range<-seq(par.aqua['d']*0.75, par.aqua['d']*1.25, length.out = sims)
  om.range<-seq(par.aqua['om']*0.75, par.aqua['om']*1.25, length.out = sims)
  gam.range<-seq(par.aqua['gam']*0.75, par.aqua['gam']*1.25, length.out = sims)
  c.range<-seq(0.015, 0.15, length.out = sims)
  p.range<-seq(5, 20, length.out = sims)
  
  paranges<-cbind(a.p=a.p.range, 
                  b.p=b.p.range, 
                  k=k.range, 
                  linf=linf.range, 
                  muP=muP.range, 
                  d=d.range, 
                  gam=gam.range, 
                  om=om.range,
                  c.range,
                  p.range)
  
#Hold all other parameters constant  
  constantparams<-matrix(ncol = length(par.aqua), nrow = sims)
  
  for(i in 1:length(par.aqua)){
    constantparams[,i] = rep(par.aqua[i],sims)
  }
  
  colnames(constantparams)<-names(par.aqua)
  vars<-colnames(paranges)

#Augment latin hypercube to sample from ###########
  #create a matrix of indices for the LHC where nrow=number of sims and ncol=number of variables
    LHC_indices<-matrix(0, nrow=sims, ncol=10)  
  #Add structure of latin hypercube
    for(j in 1:dim(LHC_indices)[2]){
      LHC_indices[,j]<-sample(1:sims, size=sims, replace=FALSE)
    }  
  #Fill LHC values
    params00<-cbind(a.p=a.p.range[LHC_indices[,1]], 
                    b.p=b.p.range[LHC_indices[,2]], 
                    k=k.range[LHC_indices[,3]], 
                    linf=linf.range[LHC_indices[,4]], 
                    muP=muP.range[LHC_indices[,5]], 
                    d=d.range[LHC_indices[,6]], 
                    gam=gam.range[LHC_indices[,7]], 
                    om=om.range[LHC_indices[,8]],
                    c=c.range[LHC_indices[,9]],
                    p=p.range[LHC_indices[,10]])
    
  params<-constantparams[, -which(colnames(constantparams) %in% vars)]
    
#Add in sampled values of parameters of interest from above
  params.fin<-cbind(params, params00)  
  
#df to fill
  params.fill = data.frame(h.t = 0,
                           h.bm = 0,
                           h.mbm = 0,
                           revenue = 0,
                           profit = 0,
                           roi = 0)
#Run model through with parameter sets from latin hypercube ######################
  for(i in 1:sims){
    parameters1<-params.fin[i,1:8]
    c = params.fin[i,9]
    p = params.fin[i,10]
    
    output = as.data.frame(ode(nstart.v, t.p, prawn_biomass, parameters1))
    output$B = ((par.aqua['a.p']*(output$L/10)^par.aqua['b.p'])/10)       
    output$Bt = output$B*output$P  
    max_Bt = max(output$B*output$P)
    
    params.fill[i,1] = output$time[output$Bt==max_Bt]   #harvest time
    params.fill[i,2] = max_Bt                       #harvest total biomass
    params.fill[i,3] = max_Bt * eta                 #harvest marketable biomass
    params.fill[i,4] = p*(params.fill[i,3]/1000)               # Raw Profit (revenue)
    params.fill[i,5] = p*(params.fill[i,3]/1000)*
                          exp(-delta*(output$time[output$Bt==max_Bt])) - 
                          c*(nstart.v["P"])  #Net profit
    params.fill[i,6] = (p*(params.fill[i,3]/1000)*
                           exp(-delta*(output$time[output$Bt==max_Bt])) - 
                           c*(nstart.v["P"])) / (c*(nstart.v["P"])) # ROI 
    if(i %% 50==0) print(i)
  }
  
#add parameters of interest to corresponding parameter sets in data frame
  params.harvest<-cbind(params.fin, params.fill)

  #save default plot settings
    opar<-par()
  
#Check scatter plots to assess sensitivity and check monotinicity ############## 
  #of harvest size and time, mean length at harvest, and pop size to tested parameters    
  par(mfrow = c(2,5), mar = c(4,3.75,1,0.4)+0.1)

  for(i in 1:length(vars)){
      plot(x = params.harvest[,i], y = params.harvest[,dim(params.harvest)[2]], pch = 16, col = 'black', cex = 0.7,
           xlab = vars[i], ylab = 'Mean length')
  }
  
  for(i in 1:length(vars)){
    plot(x = params.harvest[,i], y = params.harvest[,dim(params.harvest)[2]-1], pch = 16, col = 'green', cex = 0.7,
          xlab = vars[i], ylab = 'Pop size')
  }
  
  for(i in 1:length(vars)){
    plot(x = params.harvest[,i], y = params.harvest[,dim(params.harvest)[2]-2], pch = 16, col = 'blue', cex = 0.7,
         xlab = vars[i], ylab = 'Harvest time')
  }
  
  for(i in 1:length(vars)){
    plot(x = params.harvest[,i], y = params.harvest[,dim(params.harvest)[2]-3], pch = 16, col = 'red', cex = 0.7,
         xlab = vars[i], ylab = 'Harvest size')
  }
  
  for(i in 1:length(vars)){
    plot(x = params.harvest[,i], y = params.harvest[,dim(params.harvest)[2]-4], pch = 16, col = 'darkgreen', cex = 0.7,
         xlab = vars[i], ylab = 'Net profit')
  }
  #restore original plot settings
    par(opar)
    
#Formally analyze sensitivity of harvest time and size using PRCC ############
#net profit from harvest
  profit.pcc = pcc(X = as.data.frame(params.harvest[,c(1:length(vars))]),
                   y = as.data.frame(params.harvest$profit),
                   rank = TRUE)
    profit.pcc.df = profit.pcc$PRCC
    profit.pcc.df$var = vars
    
profit.lhsprcc = ggplot(profit.pcc.df, aes(x = var, y = original)) +
                  theme_bw()+
                  scale_y_continuous(limits = c(-0.5,0.5), breaks = c(-0.5, -0.25,0,0.25,0.5))+
                  geom_bar(fill = 'darkgreen', stat = 'identity', width = 0.25)+
                  geom_hline(yintercept = 0, lty = 2)+
                  labs(title = 'Profit sensitivity', x = 'Parameter', y = 'PRCC')
    
  profit.lhsprcc