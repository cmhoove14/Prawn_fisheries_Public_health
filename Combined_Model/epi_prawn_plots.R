#### Full snail-prawn model including epidemiological, predation, and aquaculture components

## To-do list:
  # Understand implications of diagonal transition structure
  # Simulations: basic, Gates, natural restoration (simple approximation)
  # Sensitivity analysis
  # Other modeling ideas, time permitting:
    # Incorporate seasonality?
    #incorporate finding of higher predation rate on Infected snails
#  Default settings: area = 10000, P = 5000/ha, L = 25 (0.2g post-larvae)
#  Gates settings: area = 20000; P = 500, 1000, 2500, 5000, 10000/ha; L = 67 or 100 (5g or 20g adults)

source('Combined_Model/epi_prawn_mod.R')
require(ggplot2)

## Set initial values and parameters ###########
years = 10
t.all = c(0:(years*365)+366)

#Volenhovenii stocking and harvesting events ###########
nstart.vol = c(S1 = sn.eqbm$S1, S2 = sn.eqbm$S2, S3 = sn.eqbm$S3, 
              E1 = sn.eqbm$E1, E2 = sn.eqbm$E2, E3 = sn.eqbm$E3, 
              I2 = sn.eqbm$I2, I3 = sn.eqbm$I3, W = sn.eqbm$W, 
              P = 4500, L = 25) #Optimal stocking parameters from volenhovenii ROI optimization
  harvest.t.vol = 361           #Optimal harvest time from volenhovenii optimization
  n.harvest.vol = floor(max(t.all)/harvest.t.vol)
  stocks.vol = data.frame(var = rep(c('P', 'L'), n.harvest.vol),
                          time = rep(366 + harvest.t.vol*c(0:(n.harvest.vol-1)), each = 2),
                          value = rep(c(nstart.vol['P'], nstart.vol['L']), n.harvest.vol),
                          method = rep('rep', n.harvest.vol*2))
    vol.harvests = stocks.vol[seq(1, nrow(stocks.vol), 2), 2]

  par.all.vol = c(par.aqua, par.snails,
                  a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  
#Rosenbergii stocking and harvesting events ###########
nstart.ros = c(S1 = sn.eqbm$S1, S2 = sn.eqbm$S2, S3 = sn.eqbm$S3, 
               E1 = sn.eqbm$E1, E2 = sn.eqbm$E2, E3 = sn.eqbm$E3, 
               I2 = sn.eqbm$I2, I3 = sn.eqbm$I3, W = sn.eqbm$W, 
               P = 1000, L = 25) #Optimal stocking parameters from rosenbergii ROI optimization

  harvest.t.ros = 228
  n.harvest.ros = floor(max(t.all)/harvest.t.ros)
  stocks.ros = data.frame(var = rep(c('P', 'L'), n.harvest.ros),
                          time = rep(366 + harvest.t.ros*c(0:(n.harvest.ros-1)), each = 2),
                          value = rep(c(nstart.ros['P'], nstart.ros['L']), n.harvest.ros),
                          method = rep('rep', n.harvest.ros*2))
  ros.harvests = stocks.ros[seq(1, nrow(stocks.ros), 2), 2]
  
  par.all.ros = c(par.aqua, par.snails,
                  a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
                  ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
                  th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014
  
  par.all.ros['k'] = 0.0104333333  # alternate value for M. rosenbergii, from Sampaio & Valenti 1996: 0.0104333333
  
#allow for one year run in to display epi model at eqbm ###########
t.yr1 = c(0:365) 
nstart.yr1 = c(S1 = sn.eqbm$S1, S2 = sn.eqbm$S2, S3 = sn.eqbm$S3, 
               E1 = sn.eqbm$E1, E2 = sn.eqbm$E2, E3 = sn.eqbm$E3, 
               I2 = sn.eqbm$I2, I3 = sn.eqbm$I3, W = sn.eqbm$W)

yr1 = as.data.frame(ode(nstart.yr1,t.yr1,snail_epi,par.snails))

  yr1$P = 0   #Add prawn variable to integrate with prawn stocking sims
  yr1$L = 0   #Add prawn variable to integrate with prawn stocking sims
  yr1$prev = pnbinom(2, size = 0.2, mu = yr1$W, lower.tail = FALSE)   
  yr1$S.t = (yr1$S1 + yr1$S2 + yr1$S3) / area        # density susceptible snails
  yr1$E.t = (yr1$E1 + yr1$E2 + yr1$E3) / area        # density exposed snails 
  yr1$I.t = (yr1$I2 + yr1$I3 ) / area                # density infected snails
  yr1$N.t = (yr1$S.t + yr1$E.t + yr1$I.t)            # density snails
  yr1$t.1 = (yr1$S1 + yr1$E1) / area                    # density snails of size class 1
  yr1$t.2 = (yr1$S2 + yr1$E2 + yr1$I2) / area      # density snails of size class 2
  yr1$t.3 = (yr1$S3 + yr1$E3 + yr1$I3) / area      # density snails of size class 3
  
#Simulate annual MDA for n years ########
#mda events
  mdas = data.frame(var = rep('W', years),
                    time = 365*c(1:years)+1,
                    value = rep(0.4, years),
                    method = rep('multiply', years))
  
  sim.mda = as.data.frame(ode(nstart.yr1,t.all,snail_epi,par.snails,
                          events = list(data = mdas)))
  
  sim.mda$P = 0   #Add prawn variable to integrate with prawn stocking sims
  sim.mda$L = 0   #Add prawn variable to integrate with prawn stocking sims
  sim.mda$prev = pnbinom(2, size = 0.2, mu = sim.mda$W, lower.tail = FALSE)   
  sim.mda$S.t = (sim.mda$S1 + sim.mda$S2 + sim.mda$S3) / area        # density susceptible snails
  sim.mda$E.t = (sim.mda$E1 + sim.mda$E2 + sim.mda$E3) / area        # density exposed snails 
  sim.mda$I.t = (sim.mda$I2 + sim.mda$I3 ) / area                    # density infected snails
  sim.mda$N.t = (sim.mda$S.t + sim.mda$E.t + sim.mda$I.t)            # density snails
  sim.mda$t.1 = (sim.mda$S1 + sim.mda$E1) / area                    # density snails of size class 1
  sim.mda$t.2 = (sim.mda$S2 + sim.mda$E2 + sim.mda$I2) / area      # density snails of size class 2
  sim.mda$t.3 = (sim.mda$S3 + sim.mda$E3 + sim.mda$I3) / area      # density snails of size class 3
  
  sim.mda = rbind(yr1, sim.mda)
  
#plot worm burden over time with MDA events  
  w.mda = ggplot(data = sim.mda, aes(x = time)) +
            theme_bw() +
            geom_line(aes(y = W), size = 1.25, col = 'purple') 
  
  w.mda
  
#plot snail infection dynamics over time with MDA events  
  snail.mda = ggplot(data = sim.mda, aes(x = time)) +
                theme_bw() +
                geom_line(aes(y = N.t), size = 1.25, col = 'black') +
                geom_line(aes(y = S.t), size = 1.25, col = 'green') +
                geom_line(aes(y = E.t), size = 1.25, col = 'orange') +
                geom_line(aes(y = I.t), size = 1.25, col = 'red') 
  
  snail.mda
  
#plot prevalence over time with MDA events  
  prev.mda = ggplot(data = sim.mda, aes(x = time)) +
              theme_bw() +
              geom_line(aes(y = prev), size = 1.25, col = 'red') 
  
  prev.mda
  
#plot snail size dynamics over time with MDA events  
  size.mda = ggplot(data = sim.mda, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = t.1), size = 1.25, col = 'blue') +
    geom_line(aes(y = t.2), size = 1.25, col = 'green') +
    geom_line(aes(y = t.3), size = 1.25, col = 'red') 

  size.mda
#Simulate M. volenhovenii stocking for n years ########
sim.vol = as.data.frame(ode(nstart.vol,t.all,snail_prawn_model,par.all.vol,
                              events = list(data = stocks.vol)))
  
  sim.vol$prev = pnbinom(2, size = 0.2, mu = sim.vol$W, lower.tail = FALSE)   
  sim.vol$S.t = (sim.vol$S1 + sim.vol$S2 + sim.vol$S3) / area        # density susceptible snails
  sim.vol$E.t = (sim.vol$E1 + sim.vol$E2 + sim.vol$E3) / area        # density exposed snails 
  sim.vol$I.t = (sim.vol$I2 + sim.vol$I3 ) / area                    # density infected snails
  sim.vol$N.t = (sim.vol$S.t + sim.vol$E.t + sim.vol$I.t)            # density snails
  sim.vol$t.1 = (sim.vol$S1 + sim.vol$E1) / area                          # density snails of size class 1
  sim.vol$t.2 = (sim.vol$S2 + sim.vol$E2 + sim.vol$I2) / area      # density snails of size class 2
  sim.vol$t.3 = (sim.vol$S3 + sim.vol$E3 + sim.vol$I3) / area      # density snails of size class 3
  
  sim.vol = rbind(yr1, sim.vol)
  
#plot worm burden over time with volenhovenii stocking events  
  w.vol = ggplot(data = sim.vol, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = W), size = 1.25, col = 'purple') 
  
  w.vol
  
#plot prevalence over time with volenhovenii stocking events  
  prev.vol = ggplot(data = sim.vol, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = prev), size = 1.25, col = 'red') 
  
  prev.vol
  
#plot prawn population dynamics over time with volenhovenii stocking events  
  prawn.vol = ggplot(data = sim.vol, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = P), size = 1.25, col = 'blue') 
  
  prawn.vol
  
#plot snail size dynamics over time with volenhovenii stocking events  
  size.vol = ggplot(data = sim.vol, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = t.1), size = 1.25, col = 'blue') +
    geom_line(aes(y = t.2), size = 1.25, col = 'green') +
    geom_line(aes(y = t.3), size = 1.25, col = 'red') 
  
  size.vol
#plot snail infection dynamics over time with volenhovenii stocking events  
  snail.vol = ggplot(data = sim.vol, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = N.t), size = 1.25, col = 'black') +
    geom_line(aes(y = S.t), size = 1.25, col = 'green') +
    geom_line(aes(y = E.t), size = 1.25, col = 'orange') +
    geom_line(aes(y = I.t), size = 1.25, col = 'red') 
  
  snail.vol
  
#Simulate M. rosenbergii stocking for n years ########
  sim.ros = as.data.frame(ode(nstart.ros,t.all,snail_prawn_model,par.all.ros,
                              events = list(data = stocks.ros)))
  
  sim.ros$prev = pnbinom(2, size = 0.2, mu = sim.ros$W, lower.tail = FALSE)   
  sim.ros$S.t = (sim.ros$S1 + sim.ros$S2 + sim.ros$S3) / area        # density susceptible snails
  sim.ros$E.t = (sim.ros$E1 + sim.ros$E2 + sim.ros$E3) / area        # density exposed snails 
  sim.ros$I.t = (sim.ros$I2 + sim.ros$I3 ) / area                    # density infected snails
  sim.ros$N.t = (sim.ros$S.t + sim.ros$E.t + sim.ros$I.t)            # density snails
  sim.ros$t.1 = (sim.ros$S1 + sim.ros$E1) / area                          # density snails of size class 1
  sim.ros$t.2 = (sim.ros$S2 + sim.ros$E2 + sim.ros$I2) / area      # density snails of size class 2
  sim.ros$t.3 = (sim.ros$S3 + sim.ros$E3 + sim.ros$I3) / area      # density snails of size class 3
  
  sim.ros = rbind(yr1, sim.ros)
  
#plot worm burden over time with rosenbergii stocking events  
  w.ros = ggplot(data = sim.ros, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = W), size = 1.25, col = 'purple') 
  
  w.ros
  
#plot prevalence over time with rosenbergii stocking events  
  prev.ros = ggplot(data = sim.ros, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = prev), size = 1.25, col = 'red') 
  
  prev.ros
  
#plot prawn population dynamics over time with rosenbergii stocking events  
  prawn.ros = ggplot(data = sim.ros, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = P), size = 1.25, col = 'blue') 
  
  prawn.ros
  
#plot snail infection dynamics over time with rosenbergii stocking events  
  snail.ros = ggplot(data = sim.ros, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = N.t), size = 1.25, col = 'black') +
    geom_line(aes(y = S.t), size = 1.25, col = 'green') +
    geom_line(aes(y = E.t), size = 1.25, col = 'orange') +
    geom_line(aes(y = I.t), size = 1.25, col = 'red') 
  
  snail.ros
  
#plot snail size dynamics over time with rosenbergii stocking events  
  size.ros = ggplot(data = sim.ros, aes(x = time)) +
    theme_bw() +
    geom_line(aes(y = t.1), size = 1.25, col = 'blue') +
    geom_line(aes(y = t.2), size = 1.25, col = 'green') +
    geom_line(aes(y = t.3), size = 1.25, col = 'red') 
  
  size.ros
  
  
##################
  
theme(legend.position = c(0.9, 0.5),        
                  axis.text.y = element_text(size = 12),  
                  axis.title.y = element_text(size = 15), 
                  axis.title.x = element_blank(), 
                  axis.text.x = element_blank(),
                  legend.text = element_text(size = 12),
                  legend.title = element_blank())  +    
            geom_line(aes(y = Bt/10, lty = 'Biomass'), size = 1.25) +
            labs(y = 'P', lty = "Parameter") +
            scale_x_continuous(breaks = harvests, limits = c(0, max(t.all) + 300))

harvests = stocks[seq(1, nrow(stocks), 2), 2]
harvest.labs = as.character()
  for(i in 1:n.harvest){
    harvest.labs[i] = expression(italic('T'[i]))
  }

par.all = c(par.aqua, par.snails,
            a.s = 0.187178454,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            b.s = 2.536764792,  # Allometric parameter for snail length-weight relationship, fitted to Sanna's data on B. glabrata
            ar.slp = 0.9050,    # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            #ar.int = 0.804928, # Coefficient for relationship between biomass ratio and attack rate, fitted to data from Sokolow et al. 2014
            th = 0.38561)       # Coefficient for relationship between biomass ratio and handling time, fitted to data from Sokolow et al. 2014

## Run model for n years and post process ##############
op.sp1 = as.data.frame(ode(nstart.vol, t.all, snail_prawn_model, par.all.vol,
                           events = list(data = stocks.vol)))

  op.sp1$B = ((par.all['a.p']*(op.sp1$L/10)^par.all['b.p'])/10)         # Mean prawn biomass
  op.sp1$Bt = op.sp1$B*op.sp1$P                                         # Total prawn biomass
  
#plot prawn aquaculture dynamics #########
pe.p1 = ggplot(data = op.sp1, aes(x = time)) +
          theme_bw() +
          theme(legend.position = c(0.9, 0.5),        
                axis.text.y = element_text(size = 12),  
                axis.title.y = element_text(size = 15), 
                axis.title.x = element_blank(), 
                axis.text.x = element_blank(),
                legend.text = element_text(size = 12),
                legend.title = element_blank())  +    
          geom_line(aes(y = P, lty = 'Prawns'), size = 1.25) +
          geom_line(aes(y = Bt/10, lty = 'Biomass'), size = 1.25) +
          labs(y = 'P', lty = "Parameter") +
          scale_x_continuous(breaks = harvests, limits = c(0, max(t.all) + 300))

pe.p1

#plot snail dynamics ##############
# Infection dynamics
  op.sp1$S.t = (op.sp1$S1 + op.sp1$S2 + op.sp1$S3)/area
  op.sp1$E.t = (op.sp1$E1 + op.sp1$E2 + op.sp1$E3)/area
  op.sp1$I.t = (op.sp1$I2 + op.sp1$I3)/area
  
  op.sp1$N = op.sp1$S.t + op.sp1$E.t + op.sp1$I.t
  
pe.p2 = ggplot(data = op.sp1, aes(x = time)) +
          theme_bw() +
          theme(legend.position = c(0.9, 0.5),        
                axis.text.y = element_text(size = 12),  
                axis.title.y = element_text(size = 15), 
                axis.title.x = element_blank(), 
                axis.text.x = element_blank(),
                legend.text = element_text(size = 12),
                legend.title = element_blank())  +   
          geom_line(aes(y = S.t, colour = 'S'), size = 1.25) +
          geom_line(aes(y = E.t, colour = 'E'), size = 1.25) +
          geom_line(aes(y = I.t, colour = 'I'), size = 1.25) +
          labs(y = 'Snail density', colour = "Parameter") +
          scale_color_manual(values = c('S' = 'blue', 'E' = 'orange','I' = 'red'), breaks = c('S', 'E', 'I')) +
          scale_x_continuous(breaks = harvests, limits = c(0, max(t.all) + 300))

pe.p2  
  
#plot worm burden over time ##############
pe.p3 = ggplot(data = op.sp1, aes(x = time, y = W)) +
          theme_bw() +
          theme(axis.text = element_text(size = 12),  #increase axis label size
                axis.title = element_text(size = 15), #increase axis title size
                title = element_text(size = 15)) +    #increase title size
          geom_line(color = 'purple', size = 1.25) +
          labs(x = 'time (days)', y = 'W') +
          scale_x_continuous(name = 'time (days)', 
                             breaks = c(0, harvests),
                             #labels = c(0, harvest.labs),
                             limits = c(0, max(t.all) + 300))

pe.p3  


#Combine plots to produce figure 1 #############
source('Prawn_aquaculture/prawn_aquaculture_plot.R')
      
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

fig1.layout = matrix(c(1,1,1,2,2,2,2,
                       1,1,1,2,2,2,2,
                       1,1,1,3,3,3,3,
                       4,4,4,3,3,3,3,
                       4,4,4,5,5,5,5,
                       4,4,4,5,5,5,5), ncol = 7, byrow = T)
p.blank = ggplot() + theme_bw()

windows(width = 150, height = 125)
multiplot(pa.p1, pe.p1, pe.p2, p.blank, pe.p3, layout = fig1.layout)
      
#infection dynamics as ratio of infection class to total snails ################
  op.sp1$ratioS = (op.sp1$S1 + op.sp1$S2 + op.sp1$S3)/ (op.sp1$N*area)   
  op.sp1$ratioE = (op.sp1$E1 + op.sp1$E2 + op.sp1$E3)/ (op.sp1$N*area)   
  op.sp1$ratioI = (op.sp1$I2 + op.sp1$I3)/ (op.sp1$N*area)   
  
  plot(op.sp1$time, op.sp1$ratioS, type = 'l', lwd = 2, ylim = c(0, 1), #xlim = c(0,25),
       xlab = 'time', ylab = 'snail infection class ratio',
       main = 'Snail infection dynamics')
    lines(op.sp1$time, op.sp1$ratioE, col = 'orange', lwd = 2)
    lines(op.sp1$time, op.sp1$ratioI, col = 2, lwd = 2)
    abline(v = harvests, lty = 2)
      legend('topright', legend = c('N', 'S', 'E', 'I'), col = c(1,3,'orange',2),
             lwd = 2, bty = 'n', cex = 0.8)
    
  
# Size dynamics   ##############
  op.sp1$size1 = (op.sp1$S1 + op.sp1$E1)/area
  op.sp1$size2 = (op.sp1$S2 + op.sp1$E2 + op.sp1$I2)/area
  op.sp1$size3 = (op.sp1$S3 + op.sp1$E3 + op.sp1$I3)/area
  
  plot(op.sp1$time, op.sp1$N, type = 'l', lwd = 2, ylim = c(0, max(op.sp1$N)),
       xlab = 'time', ylab = 'snail density per hectare',
       main = 'Snail size dynamics')
    lines(op.sp1$time, op.sp1$size1, col = 3, lwd = 2)
    lines(op.sp1$time, op.sp1$size2, col = 4, lwd = 2)
    lines(op.sp1$time, op.sp1$size3, col = 2, lwd = 2)
    legend('topright', legend = c('N', 'Size 1', 'Size 2', 'Size 3'), 
           col = c(1,3,4,2), lwd = 2, bty = 'n', cex = 0.8)

#snails consumed by size/infection class over time ###############
  #First get biomass ratio for each snail size class  
  ratio1 = op.sp1$B / (par.all['a.s']*0.4^par.all['b.s'])
  ratio2 = op.sp1$B / (par.all['a.s']*0.8^par.all['b.s'])
  ratio3 = op.sp1$B / (par.all['a.s']*1.2^par.all['b.s'])
  
  #get attack rates based on fit to Sanna's data, penalized by sqrt of area
  alpha1 = ifelse(ratio1 > 4.5212, par.all['ar.slp']*log(ratio1), 0)  / sqrt(area)
  alpha2 = ifelse(ratio2 > 4.5212, par.all['ar.slp']*log(ratio2), 0) / sqrt(area)
  alpha3 = ifelse(ratio3 > 4.5212, par.all['ar.slp']*log(ratio3), 0) / sqrt(area)
  
  #get handling time as sqrt of the area
  handle1 = 1/(par.all['th']*ratio1)
  handle2 = 1/(par.all['th']*ratio2)
  handle3 = 1/(par.all['th']*ratio3)
  
  #get number of snails consumed in each size/infection class by prawn population
  eatS1 = (alpha1*op.sp1$S1/area) / (1 + alpha1*handle1*op.sp1$N) * op.sp1$P
  eatS2 = (alpha2*op.sp1$S2/area) / (1 + alpha2*handle2*op.sp1$N) * op.sp1$P
  eatS3 = (alpha3*op.sp1$S3/area) / (1 + alpha3*handle3*op.sp1$N) * op.sp1$P
  
  eatE1 = (alpha1*op.sp1$E1/area) / (1 + alpha1*handle1*op.sp1$N) * op.sp1$P
  eatE2 = (alpha2*op.sp1$E2/area) / (1 + alpha2*handle2*op.sp1$N) * op.sp1$P
  eatE3 = (alpha3*op.sp1$E3/area) / (1 + alpha3*handle3*op.sp1$N) * op.sp1$P
  
  eatI2 = (alpha2*op.sp1$I2/area) / (1 + alpha2*handle2*op.sp1$N) * op.sp1$P
  eatI3 = (alpha3*op.sp1$I3/area) / (1 + alpha3*handle3*op.sp1$N) * op.sp1$P
  
  eat.total = eatS1 + eatS2 + eatS3 + eatE1 + eatE2 + eatE3 + eatI2 + eatI3
  
  #plot biomass ratio over time
  plot(op.sp1$time, ratio1, type = 'l', lwd = 2, col = 2,
       xlab = 'time', ylab = 'biomass ratio (prawn/snail)',
       ylim = c(0, max(c(ratio1, ratio2, ratio3))))
    lines(op.sp1$time, ratio2, lwd = 2, col = 3)
    lines(op.sp1$time, ratio3, lwd = 2, col = 4)
    
  #plot attack rate over time  
  plot(op.sp1$time, alpha1, type = 'l', lwd = 2, col = 2,
       xlab = 'time', ylab = 'attack rate',
       ylim = c(0, max(c(alpha1, alpha2, alpha3))))
    lines(op.sp1$time, alpha2, lwd = 2, col = 3)
    lines(op.sp1$time, alpha3, lwd = 2, col = 4) 
    
  #plot handling time over time  
  plot(op.sp1$time, handle1, type = 'l', lwd = 2, col = 2,
       xlab = 'time', ylab = 'handling time',
       ylim = c(0, max(c(handle1, handle2, handle3))))
    lines(op.sp1$time, handle2, lwd = 2, col = 3)
    lines(op.sp1$time, handle3, lwd = 2, col = 4)  
    
  #plot attack rate over biomass ratio  
  plot(ratio1, alpha1, type = 'l', lwd = 2, col = 2,
       xlab = 'biomass ratio (prawn/snail)', ylab = 'attack rate',
       ylim = c(0, max(c(alpha1, alpha2, alpha3))))
    lines(ratio2, alpha2, lwd = 2, col = 3)
    lines(ratio3, alpha3, lwd = 2, col = 4) 
    
  #plot handling time over biomass ratio  
  plot(ratio1, handle1, type = 'l', lwd = 2, col = 2,
       xlab = 'biomass ratio (prawn/snail)', ylab = 'attack rate',
       ylim = c(0, max(c(handle1, handle2, handle3))))
    lines(ratio2, handle2, lwd = 2, col = 3)
    lines(ratio3, handle3, lwd = 2, col = 4)  
    
  #plot consumption per predator over snail density
  plot(op.sp1$N, eat.total/op.sp1$P, pch = 16, cex = 0.7, 
       xlab = 'Snail density', ylab = 'Snails consumed per prawn')  
  
  #plot snail consumption of each size class across biomass ratio
  plot(ratio1, (eatS1 + eatE1)/op.sp1$P, type = 'l', lwd = 2, col = 2, ylim = c(0, max((eatS1 + eatE1)/op.sp1$P)), 
       xlab = 'biomass ratio (prawn/snail)', ylab = 'snails consumed/prawn')
    lines(ratio2, (eatS2 + eatE2 + eatI2)/op.sp1$P, lwd = 2, col = 3)
    lines(ratio3, (eatS3 + eatE3 + eatI3)/op.sp1$P, lwd = 2, col = 4)
  
  #plot consumption rate per prawn gram
  plot(op.sp1$B, eat.total/op.sp1$P, type = 'l', lwd = 2,
       xlab = 'mean prawn biomass', ylab = 'snails consumed per prawn')
    
  #plot snail consumption by predator class over time
  plot(op.sp1$time, eat.total, type = 'l', lwd = 2, ylim = c(0, max(eat.total)),
       xlab = 'time', ylab = 'snails consumed by prawn pop')
    lines(op.sp1$time, eatS1, lwd = 2, col = 3, lty = 2)
    lines(op.sp1$time, eatS2, lwd = 2, col = 3, lty = 3)
    lines(op.sp1$time, eatS3, lwd = 2, col = 3, lty = 4)
    
    lines(op.sp1$time, eatE1, lwd = 2, col = 'orange', lty = 2)
    lines(op.sp1$time, eatE2, lwd = 2, col = 'orange', lty = 3)
    lines(op.sp1$time, eatE3, lwd = 2, col = 'orange', lty = 4)
    
    lines(op.sp1$time, eatI2, lwd = 2, col = 2, lty = 3)
    lines(op.sp1$time, eatI3, lwd = 2, col = 2, lty = 4)

  #plot worm burden over time
  plot(op.sp1$time, op.sp1$W, type = 'l', col = 'purple', lwd = 2,
       xlab = 'time', ylab = 'mean worm burden')
    
    
    
    
    
    
    
# Define PZQ administration event
# Note: if administering before start of stocking, set time = 0 and instead set pzq.delay above
#   to time between MDA and first stocking
pzq.admin = 5
pzq.init = 30
pzq = data.frame(var = rep('W', pzq.admin), 
                 time = seq(pzq.init, pzq.init+(pzq.admin-1)*365, 365), 
                 value = rep(0.05, pzq.admin), 
                 method = rep('mult', pzq.admin))

# Run model and calculate outcomes of interest
output.lt = as.data.frame(ode(nstart.lt, time.lt, snail_prawn_model, parameters,
                              events = list(data = rbind(stocking, pzq))))

output.lt$S.t = output.lt$S1 + output.lt$S2 + output.lt$S3                      # Total susceptible snails
output.lt$E.t = output.lt$E1 + output.lt$E2 + output.lt$E3                      # Total exposed snails 
output.lt$I.t = output.lt$I2 + output.lt$I3                                     # Total infected snails
output.lt$N1.t = output.lt$S1 + output.lt$E1                                    # Total snails of size class 1
output.lt$N2.t = output.lt$S2 + output.lt$E2 + output.lt$I2                     # Total snails of size class 2
output.lt$N3.t = output.lt$S3 + output.lt$E3 + output.lt$I3                     # Total snails of size class 3
output.lt$N.t = output.lt$S.t + output.lt$E.t + output.lt$I.t                   # Total snails
output.lt$prev = pnbinom(2, size = 0.2, mu = output.lt$W, lower.tail = FALSE)   # Estimated prevalence, using a negative binomial dist. with k = 0.2 (fitted from EPLS data)
output.lt$B = ((parameters['a.p']*(output.lt$L/10)^parameters['b.p'])/10)       # Mean prawn biomass, transformed from length using allometric equation
output.lt$Bt = output.lt$B*output.lt$P                                          # Total prawn biomass

# Plot snail dynamics by size class
plot(x = output.lt$time, y = output.lt$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output.lt$N.t)),
     main = 'Snail Size Classes')
lines(output.lt$time, output.lt$N1.t, col = 'green', lwd = 2)
lines(output.lt$time, output.lt$N2.t, col = 'blue', lwd = 2)
lines(output.lt$time, output.lt$N3.t, col = 'red', lwd = 2)
legend('topright', legend = c('total', '1', '2', '3'), lwd = 2, col = c('black', 'green', 'blue', 'red'), cex = 0.7)

# Plot snail dynamics by infection class
plot(x = output.lt$time, y = output.lt$N.t, type = 'l', col = 'black', lwd=2, xlab = 'Time (days)', 
     ylab = 'Number of snails', ylim = c(0,max(output.lt$N.t)), 
     main = 'Snail Infection Classes')
lines(output.lt$time, output.lt$I.t, col = 'red', lwd = 2)
lines(output.lt$time, output.lt$E.t, col = 'orange', lwd = 2)
lines(output.lt$time, output.lt$S.t, col = 'green', lwd = 2)
legend('topright', legend = c('total', 'S', 'E', 'I'), lwd = 2, col = c('black', 'green', 'orange', 'red'), cex = 0.7)

# Plot mean human worm burden and prevalence of infection
plot(x = output.lt$time, y = output.lt$W, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Worm burden', ylim = c(0,max (output.lt$W)),
     main = 'Worm Burden')
plot(x = output.lt$time, y = output.lt$prev, type = 'l', col = 'red', lwd=2, xlab = 'Time (days)', 
     ylab = 'Prevalence', ylim = c(0,max (output.lt$prev)),
     main = 'Estimated Prevalence')



## Assess single-cycle outcomes over multiple stocking densities
#  NOTE: set parameters and initial values first!

# Run the model over the desired range of stocking densities (in thousands of post-larvae)
x = c(0:50)
nstart.sd = nstart
time.sd = seq(0, 500, 1)
harvest = numeric(length(x))
harvest.t = numeric(length(x))
snails = numeric(length(x))
worms = numeric(length(x))
for (i in x) {
  nstart.sd['P'] = 1000*i
  output.sd = as.data.frame(ode(nstart.sd, time.sd, snail_prawn_model, parameters,
                                events = list(data = pzq)))
  output.sd$S.t = output.sd$S1 + output.sd$S2 + output.sd$S3
  output.sd$E.t = output.sd$E1 + output.sd$E2 + output.sd$E3
  output.sd$I.t = output.sd$I2 + output.sd$I3
  output.sd$Inf.t = output.sd$E.t + output.sd$I.t
  output.sd$N.t = output.sd$S.t + output.sd$E.t + output.sd$I.t
  output.sd$B = ((parameters['a.p']*(output.sd$L/10)^parameters['b.p'])/10)
  output.sd$Bt = output.sd$B*output.sd$P
  harvest[i+1] = max(output.sd$Bt)/1000
  ht.tmp = output.sd$time[output.sd$Bt == max(output.sd$Bt)]
  harvest.t[i+1] = ifelse(ht.tmp != 0, ht.tmp, max(output.sd$time))
  snails[i+1] = output.sd$N.t[output.sd$time == harvest.t[i+1]]
  worms[i+1] = output.sd$W[output.sd$time == harvest.t[i+1]]
}

# Price estimates for profit calculation from Tamil Nadu Agricultural University, http://agritech.tnau.ac.in/fishery/fish_freshwaterprawn.html
p = 140                                           # Weighted average market price of prawns, in rupees/kg 
c = 600                                           # Cost of post-larvae, in rupees/1000 PL
delta = -log(1-0.1)/365                           # Discount rate, equivalent to 10%/year
profit = p*harvest*exp(-delta*(harvest.t)) - c*x  # Profit function, in terms of revenue (discounted by time to harvest) minus stocking costs 

# Plot estimated profit and time to harvest over stocking density
# (after one aquaculture cycle)
par(mar = c(5,5,6,5))
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     ylim = c(0, max(profit)), type = 'l', lwd = 2, main = 'Profit and time to harvest \n by stocking density')
abline(v = which.max(profit)-1, lty = 2, lwd = 2)
par(new = T)
plot(x, harvest.t, col = 'blue', axes = F, xlab = NA, ylab = NA, ylim = c(0, max(harvest.t[-1])), type = 'l', lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Harvest time (days)')
legend('topright', legend = c('Profit', 'Harvest time'), lty = 1, col = c('green', 'blue'), cex = 0.7)

# Plot estimated profit and snail abundance over stocking density
# (after one aquaculture cycle)
par(mar = c(5,5,6,5))
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     ylim = c(0, max(profit)), type = 'l', lwd = 2, main = 'Profit and snail abundance \n by stocking density')
abline(v = which.max(profit)-1, lty = 2, lwd = 2)
par(new = T)
plot(x, snails, col = 'orange', axes = F, xlab = NA, ylab = NA, ylim = c(0, max(snails)), type = 'l', lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Number of snails')
legend('topright', legend = c('Profit', 'Snails'), lty = 1, col = c('green', 'orange'), cex = 0.7)

# Plot estimated profit and mean worm burden over stocking density
# (after one aquaculture cycle)
par(mar = c(5,5,6,5))
plot(x, profit, col = 'green', xlab = 'Stocking density (thousand prawns/ha)', ylab = 'Profit (rupees/ha)',
     ylim = c(0, max(profit)), type = 'l', lwd = 2, main = 'Profit and worm burden (after PZQ) \n by stocking density')
abline(v = which.max(profit)-1, lty = 2, lwd = 2)
par(new = T)
plot(x, worms, col = 'red', axes = F, xlab = NA, ylab = NA, ylim = c(0, max(worms)), type = 'l', lwd = 2)
axis(side = 4)
mtext(side = 4, line = 3, 'Mean worm burden')
legend('topright', legend = c('Profit', 'Mean worm burden'), lty = 1, col = c('green', 'red'), cex = 0.7)



## Plots for publication

# Plot the course of the aquaculture cycle
plot.length = ggplot(output) + geom_line(aes(x=time/30, y=L/10)) + 
              theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + 
              ylab("Length (cm)")
plot.weight = ggplot(output) + geom_line(aes(x=time/30, y=B)) + 
              xlab("Time (months)") + ylab("Weight (g)")
plot.prawns = ggplot(output) + geom_line(aes(x=time/30, y=P)) + 
              theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + 
              ylab("Number of prawns")
plot.biomass = ggplot(output) + geom_line(aes(x=time/30, y=Bt/1000)) + 
               xlab("Time (months)") + ylab("Total biomass (kg)")
plot_grid(plot.length, plot.prawns, plot.weight, plot.biomass, labels = c('A', 'B', 'C', 'D'), align = 'hv')

# Plot curves for optimal harvest size and time by stocking density
plot.hsize = ggplot() + geom_line(aes(x=x, y=harvest)) + 
             xlab("Stocking density (thousand prawns/ha)") + ylab("Optimal harvest size (kg)")
plot.htime = ggplot() + geom_line(aes(x=x[-1], y=harvest.t[-1])) + 
             xlab("Stocking density (thousand prawns/ha)") + ylab("Optimal harvest time (days)") + ylim(c(0, 250))
plot_grid(plot.hsize, plot.htime, labels = c('A', 'B'), align = 'hv')

# Plot optimal profit curve
plot.profit = ggplot() + geom_line(aes(x=x, y=profit)) + geom_vline(xintercept = which.max(profit)-1, linetype = "dashed") +
              theme(axis.ticks.x = element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank()) + 
              ylab("Profit (rupees)") + ylim(c(0, max(profit)))
plot.snails = ggplot() + geom_line(aes(x=x, y=snails)) + geom_vline(xintercept = which.max(profit)-1, linetype = "dashed") +
              xlab("Stocking density (thousand prawns/ha)") + ylab("Number of snails") +
              scale_y_continuous(labels = comma)
plot_grid(plot.profit, plot.snails, labels = c('B', 'D'), align = 'hv', ncol = 1)

# Plot snail population and prevalence over 10 aquaculture cycles
plot.snails10 = ggplot(output.lt) + geom_line(aes(x=time/30, y=N.t)) + geom_line(aes(x=time/30, y=E.t+I.t), linetype = "dashed") + 
                xlab("Time (months)") + ylab("Number of snails") +
                scale_y_continuous(labels = comma)
plot.prev = ggplot(output.lt) + geom_line(aes(x=time/30, y=prev)) +
            xlab("Time (months)") + ylab("Prevalence") + ylim(c(0, max(output.lt$prev)))
plot_grid(plot.snails10, plot.prev, labels = c('A', 'B'), align = 'hv')


