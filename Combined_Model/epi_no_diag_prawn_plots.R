#### Full snail-prawn model including epidemiological, predation, and aquaculture components

## To-do list:
  # Other modeling ideas, time permitting:
    # Incorporate seasonality?
    #incorporate finding of higher predation rate on Infected snails

source('Combined_Model/epi_no_diag_prawn_mod.R')
load("Prawn_aquaculture/aquaculture_sims.RData")
require(tidyverse)

## Set initial values and parameters ###########
years = 10   #number of years to simulate
t.all = c(0:(years*365)+366)  #total time vector
cov = 0.8  #MDA coverage

#Volenhovenii stocking and harvesting events ###########
nstart.vol = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
               E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
               I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
               Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu, 
               P = opt.vol$P_nought, L = opt.vol$L_nought) #Optimal stocking parameters from volenhovenii ROI optimization

  harvest.t.vol = opt.vol$h.t           #Optimal harvest time from volenhovenii optimization
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
nstart.ros = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
               E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
               I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
               Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu, 
               P = opt.ros$P_nought, L = opt.ros$L_nought) #Optimal stocking parameters from volenhovenii ROI optimization

  harvest.t.ros = opt.ros$h.t
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
nstart.yr1 = c(S1 = allvh.eqbm$S1, S2 = allvh.eqbm$S2, S3 = allvh.eqbm$S3, 
               E1 = allvh.eqbm$E1, E2 = allvh.eqbm$E2, E3 = allvh.eqbm$E3, 
               I1 = allvh.eqbm$I1, I2 = allvh.eqbm$I2, I3 = allvh.eqbm$I3, 
               Wt = allvh.eqbm$Wt, Wu = allvh.eqbm$Wu)

yr1 = as.data.frame(ode(nstart.yr1,t.yr1,snail_epi_allvh,par.snails))

  yr1$P = 0   #Add prawn variable to integrate with prawn stocking sims
  yr1$L = 0   #Add prawn variable to integrate with prawn stocking sims
  yr1$W = cov*yr1$Wt + (1-cov)*yr1$Wu
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
  mdas = data.frame(var = rep('Wt', years),
                    time = 365*c(1:years)+1,
                    value = rep((1-cov*eff), years),
                    method = rep('multiply', years))
  
  sim.mda = as.data.frame(ode(nstart.yr1,t.all,snail_epi_allvh,par.snails,
                          events = list(data = mdas)))
  
  sim.mda$P = 0   #Add prawn variable to integrate with prawn stocking sims
  sim.mda$L = 0   #Add prawn variable to integrate with prawn stocking sims
  sim.mda$W = cov*sim.mda$Wt + (1-cov)*sim.mda$Wu
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
  w.mda = sim.mda %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    ggplot(aes(x = time, y = burden, lty = pop)) +
            theme_bw() +
            theme(axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            annotate('text', x = 0.05, y = 60, label = 'A)', size = 8) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                               labels = c(-1:years),
                               limits = c(0, (365*(years+1)+10)))
    
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
  
  sim.vol$W = cov*sim.vol$Wt + (1-cov)*sim.vol$Wu
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
  w.vol = sim.vol %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    ggplot(aes(x = time, y = burden, lty = pop)) +
              theme_bw() +
              theme(axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15)) +    #increase title size
              geom_line(size = 1.25, col = 'purple') +
              annotate('text', x = 0.05, y = 60, label = 'C)', size = 8) +
              labs(x = 'time (years)', y = expression(italic('W'))) +
              scale_y_continuous(breaks = seq(0,60,10),
                                 labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                                 limits = c(0, 61)) +
              scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                 labels = c(-1:years),
                                 limits = c(0, (365*(years+1)+10)))
  
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
snail.vol = sim.vol %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  ggplot(aes(x = time, y = density, col = class)) +
              theme_bw() +
              theme(legend.position = c(0.2, 0.75),      #place legend inside plot
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15),      #increase title size
                    legend.text = element_text(size = 14),#increase legend text size
                    legend.title = element_blank())  +    #suppress legend title
              geom_line(size = 1.25) +
              scale_color_manual(values = c('orange', 'red','black', 'green')) +
              annotate('text', x = 0.05, y = 11, label = 'E)', size = 8) +
              labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
              scale_y_continuous(breaks = seq(0,45,5),
                                 labels = seq(0,45,5),
                                 limits = c(0, 46)) +
              scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                 labels = c(-1:years),
                                 limits = c(0, (365*(years+1)+10)))
  
  snail.vol
  
#Simulate M. rosenbergii stocking for n years ########
  sim.ros = as.data.frame(ode(nstart.ros,t.all,snail_prawn_model,par.all.ros,
                              events = list(data = stocks.ros)))
  
  sim.ros$W = cov*sim.ros$Wt + (1-cov)*sim.ros$Wu
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
  w.ros = sim.ros %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    ggplot(aes(x = time, y = burden, lty = pop)) +
            theme_bw() +
            theme(axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            annotate('text', x = 0.05, y = 50, label = 'D)', size = 8) +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                               labels = c(-1:years),
                               limits = c(0, (365*(years+1)+10)))
  
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
snail.ros = sim.ros %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  ggplot(aes(x = time, y = density, col = class)) +
                    theme_bw() +
                    theme(legend.position = c(0.75, 0.75),      #place legend inside plot
                          axis.text = element_text(size = 12),  #increase axis label size
                          axis.title = element_text(size = 15), #increase axis title size
                          axis.title.x=element_blank(),         #Suppress x axis
                          axis.text.x=element_blank(),          #Suppress x axis
                          title = element_text(size = 15),      #increase title size
                          legend.text = element_text(size = 14),#increase legend text size
                          legend.title = element_blank())  +    #suppress legend title
                    geom_line(size = 1.25) +
                    annotate('text', x = 0.05, y = 11, label = 'F)', size = 8) +
                    scale_color_manual(values = c('orange', 'red','black', 'green')) +
                    labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
                    scale_y_continuous(breaks = seq(0,45,5),
                                       labels = seq(0,45,5),
                                       limits = c(0, 46)) +
                    scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                       labels = c(-1:years),
                                       limits = c(0, (365*(years+1)+10)))
  
  snail.ros
  
#Merge rosenbergii and volenhovenii stocking sims for plot of prawn pops over time ###############
  sim.ros$Species = 'M. rosenbergii'
  sim.vol$Species = 'M. volenhovenii'
  sim.rosvol = rbind(sim.vol, sim.ros)
  
p.rosvol = ggplot(sim.rosvol, aes(x = time)) +
            theme_bw() +
            theme(legend.position = c(0.8, 0.12),      #place legend inside plot
                  axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15),      #increase title size
                  legend.text = element_text(size = 14),#increase legend text size
                  legend.background = element_blank(),
                  legend.title = element_blank())  +    #suppress legend title
            geom_line(aes(y = P/area, col = Species), size = 1.25) +
            annotate('text', x = 0.05, y = 2, label = 'B)', size = 8) +
            labs(x = 'time (years)', y = expression(italic('Pm'^'-2'))) +
            scale_y_continuous(breaks = seq(0,2,0.5),
                               labels = seq(0,2,0.5),
                               limits = c(-0.2, 2.1)) +
            scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                               labels = c(-1:years),
                               limits = c(-0.2, (365*(years+1)+10)))
  p.rosvol

#Simulate annual mda along with volenhovenii stocking #############
mda.vol = rbind(mdas, stocks.vol)
  mda.vol = mda.vol[order(mda.vol$time),]
  
  sim.mda.vol = as.data.frame(ode(nstart.vol,t.all,snail_prawn_model,par.all.vol,
                              events = list(data = mda.vol)))
  
  sim.mda.vol$prev = pnbinom(2, size = 0.2, mu = sim.mda.vol$W, lower.tail = FALSE)   
  sim.mda.vol$S.t = (sim.mda.vol$S1 + sim.mda.vol$S2 + sim.mda.vol$S3) / area        # density susceptible snails
  sim.mda.vol$E.t = (sim.mda.vol$E1 + sim.mda.vol$E2 + sim.mda.vol$E3) / area        # density exposed snails 
  sim.mda.vol$I.t = (sim.mda.vol$I2 + sim.mda.vol$I3 ) / area                    # density infected snails
  sim.mda.vol$N.t = (sim.mda.vol$S.t + sim.mda.vol$E.t + sim.mda.vol$I.t)            # density snails
  sim.mda.vol$t.1 = (sim.mda.vol$S1 + sim.mda.vol$E1) / area                          # density snails of size class 1
  sim.mda.vol$t.2 = (sim.mda.vol$S2 + sim.mda.vol$E2 + sim.mda.vol$I2) / area      # density snails of size class 2
  sim.mda.vol$t.3 = (sim.mda.vol$S3 + sim.mda.vol$E3 + sim.mda.vol$I3) / area      # density snails of size class 3
  
  sim.mda.vol = rbind(yr1, sim.mda.vol)
  
#plot worm burden over time with volenhovenii stocking events and mda events
w.mda.vol = ggplot(data = sim.mda.vol, aes(x = time)) +
                    theme_bw() +
                    theme(axis.text = element_text(size = 12),  #increase axis label size
                          axis.title = element_text(size = 15), #increase axis title size
                          title = element_text(size = 15)) +    #increase title size
                    geom_line(aes(y = W), size = 1.25, col = 'purple') +
                    annotate('text', x = 0.05, y = 49, label = 'E)', size = 8) +
                    labs(x = 'time (years)', y = expression(italic('W'))) +
                    scale_y_continuous(breaks = seq(0,50,10),
                                       labels = c('0', '10', '20', '30', '40', '   50'),
                                       limits = c(0, 51)) +
                    scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                       labels = c(-1:years),
                                       limits = c(0, (365*(years+1)+10)))
  
  w.mda.vol

#Plot snail infection dynamics under repeated MDA and volenhovenii stocking
  snail.mda.vol.long = reshape(sim.mda.vol, varying = c('S.t', 'E.t', 'I.t', 'N.t'), v.names = 'infection',
                           times = c('S', 'E', 'I', 'N'), direction = 'long', drop = c('S1', 'S2', 'S3',
                                                                                       'E1', 'E2', 'E3',
                                                                                       'I2', 'I3', 'W',
                                                                                       'P', 'L', 'prev',
                                                                                       't.1', 't.2', 't.3',
                                                                                       'Species'))
  colnames(snail.mda.vol.long) = c('Infection', 'dens', 'time')
  snail.mda.vol.long$Infection = factor(snail.mda.vol.long$Infection, levels = c('N', 'S', 'E', 'I'))
  snail.mda.vol.long.sei = subset(snail.mda.vol.long, Infection != "N")
  
snail.mda.vol = ggplot(data = snail.mda.vol.long.sei, aes(x = time)) +
                  theme_bw() +
                  theme(legend.position = c(0.85, 0.55),      #place legend inside plot
                        axis.text = element_text(size = 12),  #increase axis label size
                        axis.title = element_text(size = 15), #increase axis title size
                        axis.title.x=element_blank(),         #Suppress x axis
                        axis.text.x=element_blank(),          #Suppress x axis
                        title = element_text(size = 15),      #increase title size
                        legend.text = element_text(size = 14),#increase legend text size
                        legend.background = element_blank(),  #clear legend background
                        legend.title = element_blank())  +    #suppress legend title
                  geom_line(aes(y = dens, col = Infection), size = 1.25) +
                  annotate('text', x = 0.05, y = 39, label = 'G)', size = 8) +
                  scale_color_manual(values = c('green', 'orange', 'red')) +
                  labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
                  scale_y_continuous(breaks = seq(0,40,5),
                                     labels = seq(0,40,5),
                                     limits = c(0, 40)) +
                  scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                     labels = c(-1:years),
                                     limits = c(0, (365*(years+1)+10)))
  
  snail.mda.vol
  
#Simulate annual mda along with rosenbergii stocking #############
mda.ros = rbind(mdas, stocks.ros)
  mda.ros = mda.ros[order(mda.ros$time),]
  
  sim.mda.ros = as.data.frame(ode(nstart.ros,t.all,snail_prawn_model,par.all.ros,
                                  events = list(data = mda.ros)))
  
  sim.mda.ros$prev = pnbinom(2, size = 0.2, mu = sim.mda.ros$W, lower.tail = FALSE)   
  sim.mda.ros$S.t = (sim.mda.ros$S1 + sim.mda.ros$S2 + sim.mda.ros$S3) / area        # density susceptible snails
  sim.mda.ros$E.t = (sim.mda.ros$E1 + sim.mda.ros$E2 + sim.mda.ros$E3) / area        # density exposed snails 
  sim.mda.ros$I.t = (sim.mda.ros$I2 + sim.mda.ros$I3 ) / area                    # density infected snails
  sim.mda.ros$N.t = (sim.mda.ros$S.t + sim.mda.ros$E.t + sim.mda.ros$I.t)            # density snails
  sim.mda.ros$t.1 = (sim.mda.ros$S1 + sim.mda.ros$E1) / area                          # density snails of size class 1
  sim.mda.ros$t.2 = (sim.mda.ros$S2 + sim.mda.ros$E2 + sim.mda.ros$I2) / area      # density snails of size class 2
  sim.mda.ros$t.3 = (sim.mda.ros$S3 + sim.mda.ros$E3 + sim.mda.ros$I3) / area      # density snails of size class 3
  
  sim.mda.ros = rbind(yr1, sim.mda.ros)
  
#plot worm burden over time with volenhovenii stocking events and mda events
w.mda.ros = ggplot(data = sim.mda.ros, aes(x = time)) +
                    theme_bw() +
                    theme(axis.text = element_text(size = 12),  #increase axis label size
                          axis.title = element_text(size = 15), #increase axis title size
                          title = element_text(size = 15)) +    #increase title size
                    geom_line(aes(y = W), size = 1.25, col = 'purple') +
                    annotate('text', x = 0.05, y = 49, label = 'F)', size = 8) +
                    labs(x = 'time (years)', y = expression(italic('W'))) +
                    scale_y_continuous(breaks = seq(0,50,10),
                                       labels = c('0', '10', '20', '30', '40', '   50'),
                                       limits = c(0, 51)) +
                    scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                       labels = c(-1:years),
                                       limits = c(0, (365*(years+1)+10)))
  
  w.mda.ros
  
#Plot snail infection dynamics under repeated MDA and rosenbergii stocking
snail.mda.ros.long = reshape(sim.mda.ros, varying = c('S.t', 'E.t', 'I.t', 'N.t'), v.names = 'infection',
                             times = c('S', 'E', 'I', 'N'), direction = 'long', drop = c('S1', 'S2', 'S3',
                                                                                         'E1', 'E2', 'E3',
                                                                                         'I2', 'I3', 'W',
                                                                                          'P', 'L', 'prev',
                                                                                          't.1', 't.2', 't.3',
                                                                                          'Species'))
  colnames(snail.mda.ros.long) = c('Infection', 'dens', 'time')
  snail.mda.ros.long$Infection = factor(snail.mda.ros.long$Infection, levels = c('N', 'S', 'E', 'I'))
  snail.mda.ros.long.sei = subset(snail.mda.ros.long, Infection != "N")
  
  snail.mda.ros = ggplot(data = snail.mda.ros.long.sei, aes(x = time)) +
    theme_bw() +
    theme(legend.position = "none",             #suppress legend 
          axis.text = element_text(size = 12),  #increase axis label size
          axis.title = element_text(size = 15), #increase axis title size
          axis.title.x=element_blank(),         #Suppress x axis
          axis.text.x=element_blank(),          #Suppress x axis
          title = element_text(size = 15)) +    #increase title size)      
    geom_line(aes(y = dens, col = Infection), size = 1.25) +
    annotate('text', x = 0.05, y = 39, label = 'H)', size = 8) +
    scale_color_manual(values = c('green', 'orange', 'red')) +
    labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
    scale_y_continuous(breaks = seq(0,40,5),
                       labels = seq(0,40,5),
                       limits = c(0, 40)) +
    scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                       labels = c(-1:years),
                       limits = c(0, (365*(years+1)+10)))
  
  snail.mda.ros
  
##################
##################
#plot combined figure with multiplot ###############   
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
  
  fig1.layout = matrix(c(1,2,
                         3,4,
                         5,6,
                         7,8), ncol = 2, byrow = T)

  windows(width = 1000, height = 1000)
  multiplot(w.mda, p.rosvol, 
            w.vol, w.ros, 
            w.mda.vol, w.mda.ros, 
            snail.mda.vol, snail.mda.ros, 
            layout = fig1.layout)
  
#plot alternative figure with six panels #############
  fig2.layout = matrix(c(1,2,
                         3,4,
                         5,6), ncol = 2, byrow = T)
  
  windows(width = 1000, height = 1000)
  multiplot(w.mda, p.rosvol, 
            snail.mda.vol, snail.mda.ros, 
            w.mda.vol, w.mda.ros, 
            layout = fig2.layout)
##################
##################