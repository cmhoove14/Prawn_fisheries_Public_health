#Load workspace from sims script
load("Prawn_aquaculture/aquaculture_sims.RData")

require(ggplot2)


#Plot eumetric curve for unrestricted harvest time ######
eum_dat <- rbind(opt.df, opt.df.ros)  
  
  eum_crv = ggplot(eum_dat, aes(x = P_nought/1000)) +
            theme_bw() +
            theme(legend.position = c(0.2, 0.1),        #place legend inside plot
                  axis.text = element_text(size = 15),  #increase axis label size
                  axis.title = element_text(size = 18), #increase axis title size
                  legend.text = element_text(size = 12),#increase legend text size
                  legend.title = element_blank())  +    #suppress legend title
            geom_line(aes(y = Profit, col = Species), size = 1.25) +
            #geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
            #geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
            labs(x = expression(paste('Stocking density (P', m^-2,')', sep = "")), y = 'Profit (USD)')  +
            scale_y_continuous(limits = c(-1000, 1000),
                               breaks = c(-1000, -500, 0, 500, 1000))
    eum_crv

#Plot eumetric curve for harvest time fixed at 8 months######
eum_dat8mos <- rbind(opt.df8mos, opt.df.ros8mos)  
  
  eum_crv8mos = ggplot(eum_dat8mos, aes(x = P_nought/1000)) +
            theme_bw() +
            theme(legend.position = c(0.2, 0.1),        #place legend inside plot
                  axis.text = element_text(size = 15),  #increase axis label size
                  axis.title = element_text(size = 18), #increase axis title size
                  legend.text = element_text(size = 12),#increase legend text size
                  legend.title = element_blank())  +    #suppress legend title
            geom_line(aes(y = Profit, col = Species), size = 1.25) +
            #geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
            #geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
            labs(x = expression(paste('Stocking density (P', m^-2,')', sep = "")), y = 'Profit (USD)')  +
            scale_y_continuous(limits = c(-1000, 1000),
                               breaks = c(-1000, -500, 0, 500, 1000))
    eum_crv8mos
    
#plot showing length trajectory in each species ###########
  op.spe = rbind(op.ros, op.vol)
  
    pr.l = ggplot(op.spe, aes(x = time)) +
            theme_bw() +
            theme(legend.position = c(0.2, 0.1),        #place legend inside plot
                  axis.text = element_text(size = 15),  #increase axis label size
                  axis.title = element_text(size = 18), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  legend.text = element_text(size = 12),#increase legend text size
                  legend.title = element_blank())  +    #suppress legend title
            geom_line(aes(y = L, col = Species), size = 1.25) +
            geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
            geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
            labs(x = 'time (days)', y = 'Length (mm)', col = 'Species') +
            scale_x_continuous(breaks = c(seq(0,365*2,90)),
                               labels = c(seq(0,365*2,90)),
                               limits = c(0, 365*2)) +
            scale_y_continuous(breaks = c(0, 25, 50, 100, 150,200),
                               labels = c('0', '25', '50', '100', '150', ' 200'),
                               limits = c(0, 205))  
    pr.l

#plot showing number of prawns over time ###########
    pr.P = ggplot(op.spe, aes(x = time)) +
      theme_bw() +
      theme(axis.text = element_text(size = 15),  #increase axis label size
            axis.title = element_text(size = 18), #increase axis title size
            #axis.title.x=element_blank(),         #Suppress x axis
            #axis.text.x=element_blank(),          #Suppress x axis
            #axis.ticks.x=element_blank(),        #Suppress x axis
            legend.position = 'none')  +          #suppress legend title
      geom_line(aes(y = P.dens, col = Species), size = 1.25) +
      labs(x = 'time (days)', y = expression(paste('Prawn density (Pm'^'-2',')', sep = '')), col = 'Species') +
      geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
      geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
      scale_x_continuous(breaks = c(seq(0,365*2,90)),
                         labels = c(seq(0,365*2,90)),
                         limits = c(0, 365*2)) +
      scale_y_continuous(breaks = seq(0,3, 0.5),
                         #labels = c('0', '0.5', '1.0', '1.5', '2.0'),
                         limits = c(0, 3))
    
    pr.P
    

#plot showing mean weight of each prawn over time ###########
  pr.B = ggplot(op.spe, aes(x = time)) +
      theme_bw() +
      theme(axis.text = element_text(size = 15),  #increase axis label size
            axis.title = element_text(size = 18), #increase axis title size
            axis.title.x=element_blank(),         #Suppress x axis
            axis.text.x=element_blank(),          #Suppress x axis
            legend.position = 'none')  +    #suppress legend title
      geom_line(aes(y = B, col = Species), size = 1.25) +
      geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
      geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
      labs(x = 'time (days)', y = 'Weight (g)', col = 'Species') +
      scale_x_continuous(breaks = c(seq(0,365*2,90)),
                         labels = c(seq(0,365*2,90)),
                         limits = c(0, 365*2)) +
      scale_y_continuous(breaks = seq(0,200, 50),
                         labels = c('0', '50', '100', '150', '200'),
                         limits = c(0, 200))
  pr.B
  

#plot showing total biomass of prawns over time ###########
  pr.Bt = ggplot(op.spe, aes(x = time)) +
    theme_bw() +
    theme(axis.text = element_text(size = 15),  #increase axis label size
          axis.title = element_text(size = 18), #increase axis title size
          legend.position = 'none')  +    #suppress legend title
    geom_line(aes(y = Bt, col = Species), size = 1.25) +
    geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
    geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
    labs(x = 'time (days)', y = 'Total biomass (kg)', col = 'Species') +
    scale_x_continuous(breaks = c(seq(0,365*2,90)),
                       labels = c(seq(0,365*2,90)),
                       limits = c(0, 365*2)) +
    scale_y_continuous(breaks = seq(0,200,50),
                       #labels = c('0','100','200','300','400','500'),
                       limits = c(0, 200))
    
    
  pr.Bt
 
#Combine plots to produce figure 2 #############
  
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
  
  fig2.layout = matrix(c(1,2,
                         3,4), ncol = 2, byrow = T)

  windows(width = 150, height = 100)
  multiplot(pr.l, pr.B, pr.P, pr.Bt,  layout = fig2.layout)
  
