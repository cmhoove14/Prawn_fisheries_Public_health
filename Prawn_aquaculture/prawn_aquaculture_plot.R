#Load workspace from sims script
load("Prawn_aquaculture/aquaculture_sims.RData")

require(tidyverse)
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


#Plot eumetric curve for unrestricted harvest time ######
eum_dat <- rbind(opt.df, opt.df.ros)  
  
  eum_crv = ggplot(eum_dat, aes(x = P_nought/1000)) +
              theme_bw() +
              theme(legend.position = c(0.2, 0.1),        #place legend inside plot
                    axis.text = element_text(size = 15),  #increase axis label size
                    axis.title = element_text(size = 18), #increase axis title size
                    legend.text = element_text(size = 12),#increase legend text size
                    axis.title.x = element_blank(),       #Suppress x axis since it shares it with below plot
                    axis.text.x = element_blank(),
                    legend.title = element_blank())  +    #suppress legend title
              geom_line(aes(y = Profit, col = Species), size = 1.25) +
              geom_hline(yintercept = 0, aes(col = grey20)) +
              annotate("text", x = 0, y = 1000, label = "A", size = 8) +
              labs(x = expression(paste('Stocking density (', P[0],')', sep = "")), 
                   y = expression(paste('Profit (', Pi[max], ')', sep = "")))  +
              scale_y_continuous(limits = c(-500, 1000),
                                 breaks = c(-500, 0, 500, 1000))
    eum_crv
    
eum_time = ggplot(eum_dat, aes(x = P_nought/1000)) +
            theme_bw() +
            theme(legend.position = 'none',        #place legend inside plot
                  axis.text = element_text(size = 15),  #increase axis label size
                  axis.title = element_text(size = 18), #increase axis title size
                  legend.text = element_text(size = 12),#increase legend text size
                  legend.title = element_blank())  +    #suppress legend title
            geom_line(aes(y = h.t, col = Species), size = 1.25) +
            annotate("text", x = 0, y = 730, label = "B", size = 8) +
            labs(x = expression(paste('Stocking density (', P[0],')', sep = "")), 
                 y = expression(paste('Harvest time (', T[opt]^sp, ')', sep = "")))  +
            scale_y_continuous(limits = c(0, 730),
                               breaks = c(0, 365, 730), 
                               labels = c("0", "365", "730"))  

  eum_time
  
  eum.layout = matrix(c(1,2), ncol = 1, byrow = T)

  windows(width = 100, height = 75)
  multiplot(eum_crv, eum_time, layout = eum.layout)


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
            annotate("text", x = 0, y = 200, label = "A)", size = 8) +
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
      annotate("text", x = 0, y = 3.2, label = "C)", size = 8) +
      labs(x = 'time (days)', y = expression(paste('Prawn density (Pm'^'-2',')', sep = '')), col = 'Species') +
      geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
      geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
      scale_x_continuous(breaks = c(seq(0,365*2,90)),
                         labels = c(seq(0,365*2,90)),
                         limits = c(0, 365*2)) +
      scale_y_continuous(breaks = seq(0,3, 0.5),
                         #labels = c('0', '0.5', '1.0', '1.5', '2.0'),
                         limits = c(0, 3.3))
    
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
      annotate("text", x = 0, y = 200, label = "B)", size = 8) +
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
    geom_line(aes(y = marketable_mass, col = Species), size = 1.25) +
    geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
    geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
    annotate("text", x = 0, y = 200, label = "D)", size = 8) +
    labs(x = 'time (days)', y = 'Marketable biomass (kg)', col = 'Species') +
    scale_x_continuous(breaks = c(seq(0,365*2,90)),
                       labels = c(seq(0,365*2,90)),
                       limits = c(0, 365*2)) +
    scale_y_continuous(breaks = seq(0,200,50),
                       #labels = c('0','100','200','300','400','500'),
                       limits = c(0, 210))
    
    
  pr.Bt
 
#Combine plots to produce figure 2 #############
  
  fig2.layout = matrix(c(1,2,
                         3,4), ncol = 2, byrow = T)

  windows(width = 150, height = 100)
  multiplot(pr.l, pr.B, pr.P, pr.Bt,  layout = fig2.layout)
  

  
#Produce supplementary figure of aquaculture dynamics across different stocking densities ########
  windows(width = 150, height = 150)
  aqua_sims %>% mutate(time = as.numeric(time),
                       Value = as.numeric(Value)) %>% 
    filter(Variable %in% c("P_dens", "L", "harvest_mass", "profit")) %>% 
    mutate(Variable = case_when(Variable == "P_dens" ~ "Density (P)",
                                Variable == "L" ~ "Mean Length (mm)",
                                Variable == "harvest_mass" ~ "Total biomass (kg)",
                                Variable == "profit" ~ "Profit (USD)")) %>% 
    ggplot(aes(x = time, y = Value, col = P0)) +
      geom_line(size = 1.25) +
      scale_color_manual(name = expression(P[0]),
                         values = c("gray90", "gray80", "red",
                                    "gray60", "gray50", "gray40",
                                    "gray30", "gray20")) +
      facet_grid(Variable ~ species, scales = "free_y", switch = "y") +
      theme_bw() + 
      theme(axis.text = element_text(size = 12),  #increase axis label size
            axis.title = element_text(size = 18), #increase axis title size
            legend.text = element_text(size = 15))#increase legend text size
            

