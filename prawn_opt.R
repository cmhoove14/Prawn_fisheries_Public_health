require(deSolve)

prawn_biomass=function(t, n, parameters) { 
  with(as.list(parameters),{
    
    P=n[1]
    L=n[2]
       
    Bp = (a*(L/10)^b)/10 #Mean prawn mass using conversion from length (cm) to weight (g)

    Bm = P*Bp #per prawn biomass conversion to total biomass

    dLdt= k/(1+gam*Bm)*(linf - L) #Mean prawn length growing at growth rate k limited by max length linf, and pop density (B/phi)
    
    dPdt= -P*(mu*L^-0.25 + Bm/phi) #Number of prawns subject to baseline mortality rate and density dependent mortality

    return(list(c(dPdt,dLdt)))
  }) 
} 

#Set initial values and parameters ##################
nstart=c(P=15000,L=25)
  time=seq(0,365*2,1)

#List parameters and values
parameters=c(
  a = 0.096868,
  b = 3.2944,
  gam = 1e-6, #crowding parameter that reduces growth rate at high density
  mu = 0.006136986, #baseline prawn mortality rate
  phi = 40000000, #biomass-assessed density dependence parameter in one hactare
  k = 0.00339726, #growth rate (mm/day)
  linf = 206 #max length (mm)
)

bm.start = (parameters['a']*(nstart[2]/10)^parameters['b'])/10
  bm.start

#Run & plot ############
output=as.data.frame(ode(nstart,time,prawn_biomass,parameters))
  output$B = ((parameters['a']*(output$L/10)^parameters['b'])/10)*output$P
  output$mean.size = output$B / output$P

plot(x = output$time, y = output$P, col = 'red', xlab = 'time', ylab = 'state variables', 
     type = 'l', lwd=2, xlim = c(0, max(output$time)),ylim = c(0,max(output$B/100)+100),
     main = paste('Mean start size = ', as.numeric(nstart[2]), ' mm', sep = ''))
  lines(x = output$time, y = output$B/100, col = 'blue', lwd=2)
  lines(x = output$time, y = output$mean.size, col = 'purple', lwd=2, lty=2)
  lines(x = output$time, y = output$L, col = 'green', lwd=2)
  abline(v = 8*30, lty = 2, lwd = 2)
  legend('topright', legend = c('N-prawns', 'biomass/100', 'mean size','length'), lty = 1, 
         col = c('red', 'blue', 'purple','green'), cex = 0.5)
  
start.mass = output$B[output$time==1]
  start.mass
harvest.mass = max(output$B)
  harvest.mass
harvest.size = output$mean.size[output$B==max(output$B)]  
  harvest.size
harvest.time = output$time[output$B==max(output$B)]  
  harvest.time  
  
  legend('top', legend = c(paste('starting mass=',round(start.mass), 'g', sep = ''),
                           paste('total harvest mass=',round(harvest.mass), 'g', sep = ''), 
                           paste('mean prawn mass=', round(harvest.size), 'g', sep = ''), 
                           paste('time of harvest=', round(harvest.time), 'days', sep = '')), cex=0.5)
  

mass.8mo = output$B[output$time == 30*8]
  mass.8mo
P.8mo = output$P[output$time == 30*8]
  P.8mo
L.8mo = output$L[output$time == 30*8]
  L.8mo
size.8mo = output$mean.size[output$time == 30*8]  
  size.8mo
 
  
plot(output$L/10, output$mean.size, type = 'l', xlab = 'size(cm)',
     ylab = 'mean weight (g)', xlim = c(0,25))   

plot(output$time, output$mean.size, type = 'l', xlab = 'time',
     ylab = 'mean weight (g)') 
#Fit density dependence parameters to data #####################
prawn.optim<-function(params){
  parameters['gam'] = params[1]
  parameters['phi'] = params[2]  
  
  time = seq(0,30*8,1)
  nstart1 = c(P = 5000, L = 25)
  output1=as.data.frame(ode(nstart1,time,prawn_biomass,parameters))
    mean.size1 = ((parameters['a']*(output1$L[output1$time == max(time)]/10)^parameters['b'])/10)
    surv1 = output1$P[output1$time == 0] / output1$P[output1$time == max(time)]
    
  nstart2 = c(P = 10000, L = 25)
  output2=as.data.frame(ode(nstart2,time,prawn_biomass,parameters))
    mean.size2 = ((parameters['a']*(output2$L[output2$time == max(time)]/10)^parameters['b'])/10)
    surv2 = output2$P[output2$time == 0] / output2$P[output2$time == max(time)]
    
  nstart3 = c(P = 15000, L = 25)
  output3=as.data.frame(ode(nstart3,time,prawn_biomass,parameters))
    mean.size3 = ((parameters['a']*(output3$L[output3$time == max(time)]/10)^parameters['b'])/10)
    surv3 = output3$P[output3$time == 0] / output3$P[output3$time == max(time)]
    
  nstart4 = c(P = 25000, L = 25)
  output4=as.data.frame(ode(nstart4,time,prawn_biomass,parameters))
    mean.size4 = ((parameters['a']*(output4$L[output4$time == max(time)]/10)^parameters['b'])/10)
    surv4 = output4$P[output4$time == 0] / output4$P[output4$time == max(time)]
  
  ll1 = dnorm(mean.size1, mean = 101.650, sd = 7.5067)
  ll2 = dnorm(surv1, mean = 0.69440, sd = 0.066525)
  ll3 = dnorm(mean.size2, mean = 87.730, sd = 8.1917)
  ll4 = dnorm(surv2, mean = 0.54533, sd = 0.038724)
  ll5 = dnorm(mean.size3, mean = 69.113, sd = 8.6865)
  ll6 = dnorm(surv3, mean = 0.4269, sd = 0.063707)
  ll7 = dnorm(mean.size4, mean = 55.483, sd = 3.9793)
  #ll8 = dnorm(surv4, mean = 0.28210, sd = 0.021198)

  LL.fin = ll1*ll2*ll3*ll4*ll5*ll6*ll7#*ll8
  negLL = -log(LL.fin+1)
  return(negLL)
}

prawn.optim(params = c(1e-7, 3e8))

op.df = data.frame('gams' = seq(1e-7, 1e-5, 1e-7),
                   'phis' = seq(1e7, 1e9, 1e7),
                   'negLL' = 0) 
for(i in 1:nrow(op.df)){
  params = c(op.df[i,1], op.df[i,2])
  op.df[i,3] = prawn.optim(params)
}

op.bm<-optim(par=c(1e-6, 4e8), prawn.optim, method = 'Nelder-Mead', control = list(fnscale = -1))

#Fitting to fisheries rosenbergii data from Kutty & Kurup 2010 ##############
  n.init = c(5000,10000,15000,25000)
  bm.init = ninit*0.2
  n.per.ha = c(3560,5380,5690,7020)
  mean.weight = c(101.65, 87.73,69.11, 55.48)

  plot(x = ninit, y = bm.start/10, pch = 16, col = 'blue', xlab = 'N start', ylab = 'dependent vars', ylim = c(0,520))
    points(x = ninit, y = mean.weight, pch = 16, col = 'red')
    points(x = ninit, y = survival, pch = 16, col = 'green')