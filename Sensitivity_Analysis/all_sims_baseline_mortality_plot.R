#Load workspace from sims script
rm(list = ls())

load("Prawn_aquaculture/aquaculture_sims.RData")
load("Combined_Model/epi_no_diag_prawn_sims.RData")

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

  png("~/RemaisWork/Schisto/Stanford/Prawn_Aquaculture/Figs/V3_Figs/Mortality_tests/eumetric curve/eumetric_curve_mort.png", 
      width = 750, height = 500)
  
  multiplot(eum_crv, eum_time, layout = eum.layout)

  dev.off()
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
    geom_line(aes(y = Bt, col = Species), size = 1.25) +
    geom_vline(xintercept = harvest.time.ros, lty = 2, size = 1.25) +
    geom_vline(xintercept = harvest.time.vol, lty = 3, size = 1.25) +
    annotate("text", x = 0, y = 200, label = "D)", size = 8) +
    labs(x = 'time (days)', y = 'Total biomass (kg)', col = 'Species') +
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

  png("~/RemaisWork/Schisto/Stanford/Prawn_Aquaculture/Figs/V3_Figs/Mortality_tests/aquaculture_cycle/aquaculture_summary_mort.png", 
      width = 750, height = 500)

  multiplot(pr.l, pr.B, pr.P, pr.Bt,  layout = fig2.layout)
  
  dev.off()
  
#Produce supplementary figure of aquaculture dynamics across different stocking densities ########
  png("~/RemaisWork/Schisto/Stanford/Prawn_Aquaculture/Figs/V3_Figs/Mortality_tests/stocking_density/aquaculture_sims_mort.png", 
      width = 750, height = 750)
  
  aqua_sims %>% mutate(time = as.numeric(time),
                       Value = as.numeric(Value)) %>% 
    filter(Variable %in% c("P_dens", "L", "Bt", "profit")) %>% 
    mutate(Variable = case_when(Variable == "P_dens" ~ "Density (P)",
                                Variable == "L" ~ "Mean Length (mm)",
                                Variable == "Bt" ~ "Total biomass (kg)",
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
            

  dev.off()
#plots over time with only MDA events #####
#Worm burden
  w.mda = sim.mda %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    mutate(Population = case_when(pop == "W" ~ "Mean", 
                                  pop == "Wt" ~ "Treated", 
                                  pop == "Wu" ~ "Untreated")) %>% 
    ggplot(aes(x = time, y = burden, lty = Population)) +
            theme_bw() +
            theme(legend.position = c(0.8, 0.2),
                  axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            annotate('text', x = 365, y = 60, label = 'A)', size = 8) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            geom_vline(xintercept = 365*(years+1), lty = 2)
            #scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
            #                   labels = c(-1:years),
            #                   limits = c(0, (365*(years+1)+10)))
    
  w.mda
  
#plot snail infection dynamics over time with MDA events  
snail.mda = sim.mda %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  mutate(Infection = case_when(class == "N.t" ~ "Total",
                               class == "S.t" ~ "Susceptible",
                               class == "E.t" ~ "Exposed",
                               class == "I.t" ~ "Infected")) %>% 
  mutate(Infection = factor(Infection, levels = c("Total", "Susceptible", "Exposed", "Infected"))) %>% 
  ggplot(aes(x = time, y = density, col = Infection)) +
              theme_bw() +
              theme(legend.position = c(0.75,0.4),
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15))  +   #increase title size
              geom_line(size = 1.25) +
              scale_color_manual(values = c("Total"="black", 
                                            "Susceptible"="green", 
                                            "Exposed"="orange", 
                                            "Infected"="red")) +
              annotate('text', x = 365, y = 50, label = 'B)', size = 8) +
              labs(x = 'time (years)', y = expression(paste("Snail density (", 'N'[i], 'm'^'-2', ")"))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(0, 51)) +
              geom_vline(xintercept = 365*(years+1), lty = 2)
             # scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
             #                   labels = c(-1:years),
             #                   limits = c(0, (365*(years+1)+10))) 
  
  snail.mda

#plots over time with volenhovenii stocking events ######
#worm burden
  w.vol = sim.vol %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    mutate(Population = case_when(pop == "W" ~ "Mean", 
                                  pop == "Wt" ~ "Treated", 
                                  pop == "Wu" ~ "Untreated")) %>% 
    ggplot(aes(x = time, y = burden, lty = Population)) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            annotate('text', x = 365, y = 60, label = 'C)', size = 8) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            geom_vline(xintercept = 365*(years+1), lty = 2)
            #scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
            #                   labels = c(-1:years),
            #                   limits = c(0, (365*(years+1)+10)))
  
  w.vol

#plot snail infection dynamics over time with volenhovenii stocking events  
snail.vol = sim.vol %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  mutate(Infection = case_when(class == "N.t" ~ "Total",
                               class == "S.t" ~ "Susceptible",
                               class == "E.t" ~ "Exposed",
                               class == "I.t" ~ "Infected")) %>% 
  mutate(Infection = factor(Infection, levels = c("Total", "Susceptible", "Exposed", "Infected"))) %>% 
  ggplot(aes(x = time, y = density, col = Infection)) +
              theme_bw() +
              theme(legend.position = "none",
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15))  +   #increase title size
              geom_line(size = 1.25) +
              scale_color_manual(values = c("Total"="black", 
                                            "Susceptible"="green", 
                                            "Exposed"="orange", 
                                            "Infected"="red")) +
              annotate('text', x = 365, y = 50, label = 'D)', size = 8) +
              labs(x = 'time (years)', y = expression(paste("Snail density (", 'N'[i], 'm'^'-2', ")"))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(-0.1, 51)) +
              geom_vline(xintercept = 365*(years+1), lty = 2)
             # scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
             #                   labels = c(-1:years),
             #                   limits = c(0, (365*(years+1)+10))) 
  
  snail.vol

#plots over time with rosenbergii stocking events  ###########
  w.ros = sim.ros %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    mutate(Population = case_when(pop == "W" ~ "Mean", 
                                  pop == "Wt" ~ "Treated", 
                                  pop == "Wu" ~ "Untreated")) %>% 
    ggplot(aes(x = time, y = burden, lty = Population)) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            annotate('text', x = 365, y = 60, label = 'D)', size = 8) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            geom_vline(xintercept = 365*(years+1), lty = 2)
            #scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
            #                   labels = c(-1:years),
            #                   limits = c(0, (365*(years+1)+10)))
  
  w.ros

#plot snail infection dynamics over time with rosenbergii stocking events  
snail.ros = sim.ros %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  mutate(Infection = case_when(class == "N.t" ~ "Total",
                               class == "S.t" ~ "Susceptible",
                               class == "E.t" ~ "Exposed",
                               class == "I.t" ~ "Infected")) %>% 
  mutate(Infection = factor(Infection, levels = c("Total", "Susceptible", "Exposed", "Infected"))) %>% 
  ggplot(aes(x = time, y = density, col = Infection)) +
              theme_bw() +
              theme(legend.position = "none",
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15))  +   #increase title size
              geom_line(size = 1.25) +
              scale_color_manual(values = c("Total"="black", 
                                            "Susceptible"="green", 
                                            "Exposed"="orange", 
                                            "Infected"="red")) +
              annotate('text', x = 365, y = 50, label = 'D)', size = 8) +
              labs(x = 'time (years)', y = expression(paste("Snail density (", 'N'[i], 'm'^'-2', ")"))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(-0.1, 51)) +
              geom_vline(xintercept = 365*(years+1), lty = 2)
             # scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
             #                   labels = c(-1:years),
             #                   limits = c(0, (365*(years+1)+10))) 
  
  snail.ros

#plots over time with volenhovenii stocking events and mda events ##############
w.mda.vol = sim.mda.vol %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    mutate(Population = case_when(pop == "W" ~ "Mean", 
                                  pop == "Wt" ~ "Treated", 
                                  pop == "Wu" ~ "Untreated")) %>% 
    ggplot(aes(x = time, y = burden, lty = Population)) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            annotate('text', x = 365, y = 60, label = 'E)', size = 8) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            geom_vline(xintercept = 365*(years+1), lty = 2) +
            scale_x_continuous(breaks = seq(365, 365*(years+11), 365*5),
                               labels = c(0,5,10,15,20)) 
  
  w.mda.vol

#Plot snail infection dynamics under repeated MDA and volenhovenii stocking
snail.mda.vol = sim.mda.vol %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  mutate(Infection = case_when(class == "N.t" ~ "Total",
                               class == "S.t" ~ "Susceptible",
                               class == "E.t" ~ "Exposed",
                               class == "I.t" ~ "Infected")) %>% 
  mutate(Infection = factor(Infection, levels = c("Total", "Susceptible", "Exposed", "Infected"))) %>% 
  ggplot(aes(x = time, y = density, col = Infection)) +
              theme_bw() +
              theme(legend.position = "none",
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    title = element_text(size = 15))  +   #increase title size
              geom_line(size = 1.25) +
              scale_color_manual(values = c("Total"="black", 
                                            "Susceptible"="green", 
                                            "Exposed"="orange", 
                                            "Infected"="red")) +
              annotate('text', x = 365, y = 50, label = 'F)', size = 8) +
              labs(x = 'time (years)', y = expression(paste("Snail density (", 'N'[i], 'm'^'-2', ")"))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(-0.1, 51)) +
              geom_vline(xintercept = 365*(years+1), lty = 2) +
              scale_x_continuous(breaks = seq(365, 365*(years+11), 365*5),
                                 labels = c(0, 5, 10, 15, 20))
  
  snail.mda.vol
  
#plot worm burden over time with volenhovenii stocking events and mda events
w.mda.ros = sim.mda.ros %>% dplyr::select(time, W, Wt, Wu) %>% 
    gather("pop", "burden", W:Wu) %>% 
    mutate(Population = case_when(pop == "W" ~ "Mean", 
                                  pop == "Wt" ~ "Treated", 
                                  pop == "Wu" ~ "Untreated")) %>% 
    ggplot(aes(x = time, y = burden, lty = Population)) +
            theme_bw() +
            theme(legend.position = "none",
                  axis.text = element_text(size = 12),  #increase axis label size
                  axis.title = element_text(size = 15), #increase axis title size
                  axis.title.x=element_blank(),         #Suppress x axis
                  axis.text.x=element_blank(),          #Suppress x axis
                  title = element_text(size = 15)) +    #increase title size
            geom_line(size = 1.25, col = 'purple') +
            labs(x = 'time (years)', y = expression(italic('W'))) +
            annotate('text', x = 365, y = 60, label = 'E)', size = 8) +
            scale_y_continuous(breaks = seq(0,60,10),
                               labels = c('0', '10', '20', '30', '40', '   50', '   60'),
                               limits = c(0, 61)) +
            geom_vline(xintercept = 365*(years+1), lty = 2)
            #scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
            #                   labels = c(-1:years),
            #                   limits = c(0, (365*(years+1)+10)))
  
  w.mda.ros
  
#Plot snail infection dynamics under repeated MDA and rosenbergii stocking
snail.mda.ros = sim.mda.ros %>% dplyr::select(time, S.t, E.t, I.t, N.t) %>% 
  gather("class", "density", S.t:N.t) %>% 
  mutate(Infection = case_when(class == "N.t" ~ "Total",
                               class == "S.t" ~ "Susceptible",
                               class == "E.t" ~ "Exposed",
                               class == "I.t" ~ "Infected")) %>% 
  mutate(Infection = factor(Infection, levels = c("Total", "Susceptible", "Exposed", "Infected"))) %>% 
  ggplot(aes(x = time, y = density, col = Infection)) +
              theme_bw() +
              theme(legend.position = "none",
                    axis.text = element_text(size = 12),  #increase axis label size
                    axis.title = element_text(size = 15), #increase axis title size
                    axis.title.x=element_blank(),         #Suppress x axis
                    axis.text.x=element_blank(),          #Suppress x axis
                    title = element_text(size = 15))  +   #increase title size
              geom_line(size = 1.25) +
              scale_color_manual(values = c("Total"="black", 
                                            "Susceptible"="green", 
                                            "Exposed"="orange", 
                                            "Infected"="red")) +
              annotate('text', x = 365, y = 50, label = 'F)', size = 8) +
              labs(x = 'time (years)', y = expression(paste("Snail density (", 'N'[i], 'm'^'-2', ")"))) +
              scale_y_continuous(breaks = seq(0,50,10),
                                 labels = seq(0,50,10),
                                 limits = c(-0.1, 51)) +
              geom_vline(xintercept = 365*(years+1), lty = 2)
             # scale_x_continuous(breaks = seq(0, 365*(years+1), 365),
             #                   labels = c(-1:years),
             #                   limits = c(0, (365*(years+1)+10))) 
  
  snail.mda.ros

#Plot multipanel epi figure #######
  
  fig1.layout = matrix(c(1,2,
                         3,4,
                         5,6), ncol = 2, byrow = T)

  png("~/RemaisWork/Schisto/Stanford/Prawn_Aquaculture/Figs/V3_Figs/Mortality_tests/epi_sims/epi_sims_prawn_mort.png", 
      width = 1000, height = 750)
  
  multiplot(w.mda, snail.mda, 
            w.ros, snail.ros, 
            w.mda.ros, snail.mda.ros,
            layout = fig1.layout)
  
  dev.off()
            