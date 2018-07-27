load("Combined_Model/epi_no_diag_prawn_sims.RData")
require(tidyverse)

#plots over time with only MDA events #####
#Worm burden
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
snail.mda = sim.mda %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  ggplot(aes(x = time, y = density, col = class)) +
              theme_bw() +
              theme(legend.position = c(0.9, 0.9),      #place legend inside plot
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15),      #increase title size
                    legend.text = element_text(size = 14),#increase legend text size
                    legend.title = element_blank())  +    #suppress legend title
              geom_line(size = 1.25) +
              scale_color_manual(values = c('orange', 'red','black', 'green')) +
              annotate('text', x = 0.05, y = 50, label = 'E)', size = 8) +
              labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(0, 51)) +
              scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                 labels = c(-1:years),
                                 limits = c(0, (365*(years+1)+10))) 
  
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
  
  
#plots over time with volenhovenii stocking events ######
#worm burden
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
              theme(legend.position = c(0.9, 0.9),      #place legend inside plot
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15),      #increase title size
                    legend.text = element_text(size = 14),#increase legend text size
                    legend.title = element_blank())  +    #suppress legend title
              geom_line(size = 1.25) +
              scale_color_manual(values = c('orange', 'red','black', 'green')) +
              annotate('text', x = 0.05, y = 50, label = 'E)', size = 8) +
              labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(-.10, 51)) +
              scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                 labels = c(-1:years),
                                 limits = c(0, (365*(years+1)+10)))
  
  snail.vol
  
  
#plots over time with rosenbergii stocking events  ###########
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
                    theme(legend.position = c(0.9, 0.9),      #place legend inside plot
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
                    scale_y_continuous(breaks = seq(0,50,10),
                                       labels = seq(0,50,10),
                                       limits = c(-.10, 51)) +
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
            scale_y_continuous(breaks = seq(0,3,0.5),
                               labels = seq(0,3,0.5),
                               limits = c(-0.2, 3.1)) +
            scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                               labels = c(-1:years),
                               limits = c(-0.2, (365*(years+1)+10)))
  p.rosvol

  
#plots over time with volenhovenii stocking events and mda events ##############
w.mda.vol = sim.mda.vol %>% dplyr::select(time, W, Wt, Wu) %>% 
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
  
  w.mda.vol

#Plot snail infection dynamics under repeated MDA and volenhovenii stocking
snail.mda.vol = sim.mda.vol %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  ggplot(aes(x = time, y = density, col = class)) +
                  theme_bw() +
                  theme(legend.position = c(0.9, 0.9),      #place legend inside plot
                        axis.text = element_text(size = 12),  #increase axis label size
                        axis.title = element_text(size = 15), #increase axis title size
                        axis.title.x=element_blank(),         #Suppress x axis
                        axis.text.x=element_blank(),          #Suppress x axis
                        title = element_text(size = 15),      #increase title size
                        legend.text = element_text(size = 14),#increase legend text size
                        legend.background = element_blank(),  #clear legend background
                        legend.title = element_blank())  +    #suppress legend title
                  geom_line(size = 1.25) +
                  annotate('text', x = 0.05, y = 50, label = 'G)', size = 8) +
                  scale_color_manual(values = c('orange', 'red','black', 'green')) +
                  labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
                  scale_y_continuous(breaks = seq(0,50,10),
                                     labels = seq(0,50,10),
                                     limits = c(0, 51)) +
                  scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
                                     labels = c(-1:years),
                                     limits = c(0, (365*(years+1)+10)))
  
  snail.mda.vol
  
  
#plot worm burden over time with volenhovenii stocking events and mda events
w.mda.ros = sim.mda.ros %>% dplyr::select(time, W, Wt, Wu) %>% 
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
  
  w.mda.ros
  
#Plot snail infection dynamics under repeated MDA and rosenbergii stocking
snail.mda.ros = sim.mda.ros %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  ggplot(aes(x = time, y = density, col = class)) +
    theme_bw() +
    theme(legend.position = "none",             #suppress legend 
          axis.text = element_text(size = 12),  #increase axis label size
          axis.title = element_text(size = 15), #increase axis title size
          axis.title.x=element_blank(),         #Suppress x axis
          axis.text.x=element_blank(),          #Suppress x axis
          title = element_text(size = 15)) +    #increase title size)      
    geom_line(size = 1.25) +
    annotate('text', x = 0.05, y = 50, label = 'H)', size = 8) +
    scale_color_manual(values = c('orange', 'red','black', 'green')) +
    labs(x = 'time (years)', y = expression(paste('N'[i], 'm'^'-2'))) +
    scale_y_continuous(breaks = seq(0,50,10),
                       labels = seq(0,50,10),
                       limits = c(0, 51)) +
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