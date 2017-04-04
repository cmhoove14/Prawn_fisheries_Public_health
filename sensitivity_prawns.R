source('snail_prawn_full.R')
require(ggplot2)
require(sensitivity)
require(Rmisc)

#Start with parameters that affect prawn aquaculture dynamics ##################
  #Note: these parameter ranges are completely arbitrary;
  #should be updated with more informed max/min values and plausible distributions
  sims=250

#Prawn parameters
  a.p.range<-seq(parameters['a.p']-0.03, parameters['a.p']+0.03, length.out = sims)
  b.p.range<-seq(parameters['b.p']-0.5, parameters['b.p']+0.5, length.out = sims)
  k.range<-seq(parameters['k']/10, parameters['k']*10, length.out = sims)
  linf.range<-seq(parameters['linf']-30, parameters['linf']+30, length.out = sims)
  muP.range<-seq(parameters['muP']/10, parameters['muP']*3, length.out = sims)
  d.range<-seq(parameters['d']/10, parameters['d']*10, length.out = sims)
  phi.range<-seq(parameters['phi']/100, parameters['phi']*100, length.out = sims)
  gam.range<-seq(parameters['gam']/100, parameters['gam']*100, length.out = sims)
  
  paranges<-cbind(a.p=a.p.range, 
                  b.p=b.p.range, 
                  k=k.range, 
                  linf=linf.range, 
                  muP=muP.range, 
                  d=d.range, 
                  gam=gam.range, 
                  phi=phi.range)
  
#Hold all other parameters constant  
  constantparams<-matrix(ncol = length(parameters), nrow = sims)
  
  for(i in 1:length(parameters)){
    constantparams[,i] = rep(parameters[i],sims)
  }
  
  colnames(constantparams)<-names(parameters)
  vars<-colnames(paranges)

#First check for monotinicity between variables and outcomes #########
  outputfill1<-matrix(ncol = length(vars), nrow = sims)
  outputfill2<-matrix(ncol = length(vars), nrow = sims)
  outputfill3<-matrix(ncol = length(vars), nrow = sims)
  outputfill4<-matrix(ncol = length(vars), nrow = sims)
  
  
  for(j in 1:length(vars)){
    for(i in 1:sims){
      print(c(j, i))
      
      other<-constantparams[, -which(colnames(constantparams) %in% vars[j])]
      test<-paranges[, which(colnames(paranges) %in% vars[j])]
      
      parametersuse<-cbind(other, test) 
      colnames(parametersuse)[dim(parametersuse)[2]]<-vars[j]
      
      parameters1<-parametersuse[i,]
      
      output = as.data.frame(ode(nstart, time, snail_prawn_model, parameters1))
      output$B = ((parameters['a.p']*(output$L/10)^parameters['b.p'])/10)       
      output$Bt = output$B*output$P                                             
      harvest.size = output$B[output$Bt == max(output$Bt)]                      
      harvest.time = output$time[output$Bt == max(output$Bt)]    
      outputfill1[i,j] = harvest.size
      outputfill2[i,j] = harvest.time
      outputfill3[i,j] = output$P[output$Bt == max(output$Bt)] 
      outputfill4[i,j] = output$L[output$Bt == max(output$Bt)]
      print(c(i, j, outputfill1[i,j], outputfill2[i,j],  outputfill3[i,j], outputfill4[i,j]))
      
      
    }
    par(mfrow = c(4,1), mar = c(4,3.75,1,0.4)+0.1)
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill1[,j], 
         xlab = vars[j], ylab = 'harvest size',
         pch = 16, cex = 0.75, col = 'red', 
         ylim = c(0,120))
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill2[,j], 
         xlab = vars[j], ylab = 'harvest time',
         pch = 16, cex = 0.75, col = 'blue', 
         ylim = c(0,365*2))
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill3[,j], 
         xlab = vars[j], ylab = 'Population at harvest',
         pch = 16, cex = 0.75, col = 'green', 
         ylim = c(0,5000))
    
    plot(x = parametersuse[,dim(parametersuse)[2]], y = outputfill4[,j], 
         xlab = vars[j], ylab = 'Mean size at harvest',
         pch = 16, cex = 0.75, col = 'grey50', 
         ylim = c(0,150))
  }
  
#Change parameter ranges based on results from above ###########
  #make k range from 0.02 to previous max
    k.range<-seq(0.02, max(k.range), length.out = sims)
  #make d range from -3 to -0.03
    d.range<-seq(-3, -0.03, length.out = sims)

#Augment latin hypercube to sample from ###########
  #create a matrix of indices for the LHC where nrow=number of sims and ncol=number of variables
    LHC_indices<-matrix(0, nrow=sims, ncol=8)  
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
                    phi=phi.range[LHC_indices[,8]])
    
  params<-constantparams[, -which(colnames(constantparams) %in% vars)]
    
  #Add in sampled values of parameters of interest from above
    params.fin<-cbind(params, params00)  
#Run model through with parameter sets from latin hypercube ######################
  hs<-numeric()
  ht<-numeric()
  pp<-numeric()
  ll<-numeric()
  
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
    pp[i] = output$P[output$Bt == max(output$Bt)] 
    ll[i] = output$L[output$Bt == max(output$Bt)]
    print(c(i, hs[i], ht[i], pp[i], ll[i]))
  }
  
#add parameters of interest to corresponding parameter sets in data frame
  params.harvest<-as.data.frame(cbind(params.fin, hs, ht, pp, ll))
  ph.test<-cbind(params.harvest[, -which(colnames(params.harvest) %in% colnames(params))])
  
  #save default plot settings
    opar<-par()
  
#Check scatter plots to assess sensitivity ############## 
  #of harvest size and time, mean length at harvest, and pop size to tested parameters    
  par(mfrow = c(2,4), mar = c(4,3.75,1,0.4)+0.1)

  for(i in 1:length(vars)){
      plot(x = ph.test[,i], y = ph.test[,dim(ph.test)[2]], pch = 16, col = 'black', cex = 0.7,
           xlab = vars[i], ylab = 'Mean length')
  }
  
  for(i in 1:length(vars)){
    plot(x = ph.test[,i], y = ph.test[,dim(ph.test)[2]-1], pch = 16, col = 'green', cex = 0.7,
          xlab = vars[i], ylab = 'Pop size')
  }
  
  for(i in 1:length(vars)){
    plot(x = ph.test[,i], y = ph.test[,dim(ph.test)[2]-2], pch = 16, col = 'blue', cex = 0.7,
         xlab = vars[i], ylab = 'Harvest time')
  }
  
  for(i in 1:length(vars)){
    plot(x = ph.test[,i], y = ph.test[,dim(ph.test)[2]-3], pch = 16, col = 'red', cex = 0.7,
         xlab = vars[i], ylab = 'Harvest size')
  }
  
  #restore original plot settings
    par(opar)
    
#Formally analyze sensitivity of harvest time and size using PRCC ############
#harvest time prcc    
  ht.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
              y = as.data.frame(ph.test[,c(dim(ph.test)[2]-2)]),
              rank = TRUE)
  
  ht.pcc.df<-ht.pcc$PRCC
  ht.pcc.df$var = vars

  htp = ggplot(ht.pcc.df, aes(x = var, y = original)) +
          theme_bw()+
          scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
          geom_bar(fill = 'blue', stat = 'identity', width = 0.25)+
          geom_hline(yintercept = 0, lty = 2)+
          labs(title = 'Harvest time sensitivity', x = 'Parameter', y = 'PRCC')
  
#harvest size prcc    
  hs.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
              y = as.data.frame(ph.test[,c(dim(ph.test)[2]-3)]),
              rank = TRUE)
  
  hs.pcc.df<-hs.pcc$PRCC
  hs.pcc.df$var = vars

  hsp = ggplot(hs.pcc.df, aes(x = var, y = original)) +
          theme_bw()+
          scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
          geom_bar(fill = 'red', stat = 'identity', width = 0.25)+
          geom_hline(yintercept = 0, lty = 2)+
          labs(title = 'Harvest size sensitivity', x = 'Parameter', y = 'PRCC')

#prawn population prcc
  pp.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
              y = as.data.frame(ph.test[,c(dim(ph.test)[2]-1)]),
              rank = TRUE)
  
  pp.pcc.df<-pp.pcc$PRCC
  pp.pcc.df$var = vars

  ppp = ggplot(pp.pcc.df, aes(x = var, y = original)) +
          theme_bw()+
          scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
          geom_bar(fill = 'green', stat = 'identity', width = 0.25)+
          geom_hline(yintercept = 0, lty = 2)+
          labs(title = 'Prawn population sensitivity', x = 'Parameter', y = 'PRCC')

#mean length prcc
  ll.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
              y = as.data.frame(ph.test[,c(dim(ph.test)[2])]),
              rank = TRUE)
  
  ll.pcc.df<-ll.pcc$PRCC
  ll.pcc.df$var = vars

  llp = ggplot(ll.pcc.df, aes(x = var, y = original)) +
          theme_bw()+
          scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
          geom_bar(fill = 'black', stat = 'identity', width = 0.25)+
          geom_hline(yintercept = 0, lty = 2)+
          labs(title = 'Mean length sensitivity', x = 'Parameter', y = 'PRCC')
  
  multiplot(htp, hsp, ppp, llp, layout = matrix(c(1,2,3,4), nrow = 2, byrow = T))
  
#See how stocking size affects sensitivity ##########  
  hs<-numeric()
  ht<-numeric()
  pp<-numeric()
  ll<-numeric()
  
  dens.range<-seq(1000, 20000, by = 500)
  
  hs.fill<-matrix(nrow = ncol(LHC_indices), ncol = length(dens.range))
  ht.fill<-matrix(nrow = ncol(LHC_indices), ncol = length(dens.range))
  pp.fill<-matrix(nrow = ncol(LHC_indices), ncol = length(dens.range))
  ll.fill<-matrix(nrow = ncol(LHC_indices), ncol = length(dens.range))
  
  area = 10000
  time1 = seq(0, 365*2, 1)
  
  for(j in 1:length(dens.range)){
    for(i in 1:sims){
    parameters1<-params.fin[i,]
    nstart1= c(S1 = 3.54*area, S2 = 1.02*area, S3 = 0.22*area, 
               E1 = 3.51*area, E2 = 1.1*area, E3 = 0.25*area, 
               I2 = 1.13*area, I3 = 0.23*area, 
               W = 72, 
               P = dens.range[j]*area/10000, L = 25)
    
    output = as.data.frame(ode(nstart1, time1, snail_prawn_model, parameters1))
    output$B = ((parameters['a.p']*(output$L/10)^parameters['b.p'])/10)       
    output$Bt = output$B*output$P   
    
    harvest.size = output$B[output$Bt == max(output$Bt)]                      
    harvest.time = output$time[output$Bt == max(output$Bt)]
    
    hs[i] = harvest.size
    ht[i] = harvest.time
    pp[i] = output$P[output$Bt == max(output$Bt)] 
    ll[i] = output$L[output$Bt == max(output$Bt)]
    print(c(j, i, hs[i], ht[i], pp[i], ll[i]))
    
  #add parameters of interest to corresponding parameter sets in data frame
    params.harvest<-as.data.frame(cbind(params.fin, hs, ht, pp, ll))
    ph.test<-cbind(params.harvest[, -which(colnames(params.harvest) %in% colnames(params))])
    
    #harvest size
    hs.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
                y = as.data.frame(ph.test[,c(dim(ph.test)[2]-3)]),
                rank = TRUE)
    
    hs.fill[,j]<-hs.pcc$PRCC[[1]]
    
    #harvest time
    ht.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
                y = as.data.frame(ph.test[,c(dim(ph.test)[2]-2)]),
                rank = TRUE)
    
    ht.fill[,j]<-ht.pcc$PRCC[[1]]
    
    #harvest population
    pp.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
                y = as.data.frame(ph.test[,c(dim(ph.test)[2]-1)]),
                rank = TRUE)
    
    pp.fill[,j]<-pp.pcc$PRCC[[1]]
    
    #harvest mean prawn size
    ll.pcc<-pcc(X = as.data.frame(ph.test[,c(1:length(vars))]),
                y = as.data.frame(ph.test[,c(dim(ph.test)[2])]),
                rank = TRUE)
    
    ll.fill[,j]<-ll.pcc$PRCC[[1]]
  }
  
}
#Post process ############  
  hsdf<-as.data.frame(hs.fill)  
    colnames(hsdf)<-dens.range
    hsdf$mean<-numeric(length=length(vars))
    hsdf$sd<-numeric(length=length(vars))
  
  for(i in 1:nrow(hsdf)){
    hsdf[i, length(dens.range)+1] = mean(as.numeric(hsdf[i, c(1:length(dens.range))]))
    hsdf[i, length(dens.range)+2] = sd(as.numeric(hsdf[i, c(1:length(dens.range))]))
  }
    hsdf$vars<-vars
    
  htdf<-as.data.frame(ht.fill)  
    colnames(htdf)<-dens.range
    htdf$mean<-numeric(length=length(vars))
    htdf$sd<-numeric(length=length(vars))
    
  for(i in 1:nrow(htdf)){
    htdf[i, length(dens.range)+1] = mean(as.numeric(htdf[i, c(1:length(dens.range))]))
    htdf[i, length(dens.range)+2] = sd(as.numeric(htdf[i, c(1:length(dens.range))]))
  }
    htdf$vars<-vars
  
  ppdf<-as.data.frame(pp.fill)  
    colnames(ppdf)<-dens.range
    ppdf$mean<-numeric(length=length(vars))
    ppdf$sd<-numeric(length=length(vars))
    
  for(i in 1:nrow(ppdf)){
    ppdf[i, length(dens.range)+1] = mean(as.numeric(ppdf[i, c(1:length(dens.range))]))
    ppdf[i, length(dens.range)+2] = sd(as.numeric(ppdf[i, c(1:length(dens.range))]))
  }
    ppdf$vars<-vars
  
  lldf<-as.data.frame(ll.fill)  
    colnames(lldf)<-dens.range
    lldf$mean<-numeric(length=length(vars))
    lldf$sd<-numeric(length=length(vars))
    
  for(i in 1:nrow(lldf)){
    lldf[i, length(dens.range)+1] = mean(as.numeric(lldf[i, c(1:length(dens.range))]))
    lldf[i, length(dens.range)+2] = sd(as.numeric(lldf[i, c(1:length(dens.range))]))
  }
    lldf$vars<-vars
#Static mean +/- sd of prcc plots ############
    ggplot(hsdf, aes(x = vars, y = mean)) +
      theme_bw()+
      scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
      geom_bar(fill = 'red', stat = 'identity', width = 0.25)+
      labs(title = 'Harvest size sensitivity', x = 'Parameter', y = 'Mean +/- sd PRCC')+
      geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.15)
    
    ggplot(htdf, aes(x = vars, y = mean)) +
      theme_bw()+
      scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
      geom_bar(fill = 'blue', stat = 'identity', width = 0.25)+
      labs(title = 'Harvest time sensitivity', x = 'Parameter', y = 'Mean +/- sd PRCC')+
      geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.15)
    
    ggplot(ppdf, aes(x = vars, y = mean)) +
      theme_bw()+
      scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
      geom_bar(fill = 'green', stat = 'identity', width = 0.25)+
      labs(title = 'Prawn pop sensitivity', x = 'Parameter', y = 'Mean +/- sd PRCC')+
      geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.15)
    
    ggplot(lldf, aes(x = vars, y = mean)) +
      theme_bw()+
      scale_y_continuous(limits = c(-1,1), breaks = c(-1,0,1))+
      geom_bar(fill = 'grey50', stat = 'identity', width = 0.25)+
      labs(title = 'Mean length sensitivity', x = 'Parameter', y = 'Mean +/- sd PRCC')+
      geom_errorbar(aes(x = vars, ymin = (mean - sd), ymax = (mean + sd)), width = 0.15)
 
  
    
#See how sensitivity varies based on starting density ############
  plot(x = as.numeric(colnames(hsdf)[c(1:length(dens.range))]), 
       y = hsdf[1,c(1:length(dens.range))], type = 'l', lwd = 2,
       ylim = c(-1,1), xlab = 'Starting density', ylab = 'prcc',
       main = 'Harvest size sensitivity over starting density')
    for(i in 2:nrow(hsdf)){
      lines(x = as.numeric(colnames(hsdf)[1:length(dens.range)]), 
            y = hsdf[i,c(1:length(dens.range))], col = i, lwd = 2)
    }
    legend('bottomright', legend = vars, col = c(1:8), lwd=1, cex=0.5)
    
  plot(x = as.numeric(colnames(htdf)[c(1:length(dens.range))]), 
       y = htdf[1,c(1:length(dens.range))], type = 'l', lwd = 2,
       ylim = c(-1,1), xlab = 'Starting density', ylab = 'prcc',
       main = 'Harvest time sensitivity over starting density')
    for(i in 2:nrow(htdf)){
      lines(x = as.numeric(colnames(htdf)[1:length(dens.range)]), 
            y = htdf[i,c(1:length(dens.range))], col = i, lwd = 2)
    }
    legend('bottomright', legend = vars, col = c(1:8), lwd=1, cex=0.5) 
    
  plot(x = as.numeric(colnames(ppdf)[c(1:length(dens.range))]), 
       y = ppdf[1,c(1:length(dens.range))], type = 'l', lwd = 2,
       ylim = c(-1,1), xlab = 'Starting density', ylab = 'prcc',
       main = 'Prawn pop at harvest sensitivity over starting density')
    for(i in 2:nrow(ppdf)){
      lines(x = as.numeric(colnames(ppdf)[1:length(dens.range)]), 
            y = ppdf[i,c(1:length(dens.range))], col = i, lwd = 2)
    }
    legend('bottomright', legend = vars, col = c(1:8), lwd=1, cex=0.5) 
    
  plot(x = as.numeric(colnames(lldf)[c(1:length(dens.range))]), 
       y = lldf[1,c(1:length(dens.range))], type = 'l', lwd = 2,
       ylim = c(-1,1), xlab = 'Starting density', ylab = 'prcc',
       main = 'Mean length at harvest sensitivity over starting density')
    for(i in 2:nrow(lldf)){
      lines(x = as.numeric(colnames(lldf)[1:length(dens.range)]), 
            y = htdf[i,c(1:length(dens.range))], col = i, lwd = 2)
    }
    legend('bottomright', legend = vars, col = c(1:8), lwd=1, cex=0.5) 
    
    