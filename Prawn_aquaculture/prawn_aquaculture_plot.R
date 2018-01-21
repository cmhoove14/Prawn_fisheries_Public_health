#Generate plot of aquaculture cycle over two years to display in model summary figure
source('Prawn_aquaculture/prawn_aquaculture_mod.R')
source('Prawn_aquaculture/macrobrachium_aquaculture_data.R')

t.p = c(0:(365*2))
#run model through 2 years for M. volenhovenii ######
par.aqua['k'] = 0.00339726  # alternate value for M. rosenbergii, from Sampaio & Valenti 1996
nstart.p.vol = c(P = 8000, L = 38)   #optimal starting conditions for M. volenhovenii as estimated in prawn_aquaculture_opt.R script
op.vol = as.data.frame(ode(nstart.p.vol,t.p,prawn_biomass,par.aqua))
  c = 0.1+0.005*(nstart.p.vol["L"]-38) # make juvenile prawns cost a function of their size with reference at $0.1/juvenile

#post-process to estimate additional parameters
# 0) prawn density rather than raw prawn numbers
  op.vol$P.dens = op.vol$P / 10000
# 1) mean prawn biomass (allometric function)
  op.vol$B = ((par.aqua['a.p']*(op.vol$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op.vol$Bt = op.vol$B*op.vol$P / 1000
# 2.5) marketable prawns at harvest
  eta.vol = predict(eta.lm, newdata = data.frame(dens = nstart.p.vol["P"]/area)) # Fraction of harvest that's marketable as function of stocking density
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op.vol$profit = eta.vol*p*(op.vol$Bt/1000)*exp(-delta*(op.vol$t)) - c*(nstart.p["P"])
# 4) starting total biomass
  start.mass.kg.vol = op.vol$Bt[op.vol$time==0]/1000
# 5) harvest mass in kg (harvest assumed to occur when biomass is maximized)   
  harvest.mass.kg.vol = eta.vol*max(op.vol$Bt)/1000
# 6) average mass of prawns at harvest  
  harvest.size.vol = op.vol$B[op.vol$Bt==max(op.vol$Bt)]
# 7) time of harvest (when biomass is maximized)
  harvest.time.vol = op.vol$time[op.vol$Bt==max(op.vol$Bt)]
# 8) time in months
  op.vol$time.mo = op.vol$time / 30
# 9) species
  op.vol$Species = 'M. volenhovenii'

#run model through 2 years for M. rosenbergii ######
par.aqua['k'] = 0.0104333333  # alternate value for M. rosenbergii, from Sampaio & Valenti 1996
nstart.p.ros = c(P = 19000, L = 38)   #optimal starting conditions for M. rosenbergii as estimated in prawn_aquaculture_opt.R script
op.ros = as.data.frame(ode(nstart.p.ros,t.p,prawn_biomass,par.aqua))
  c = 0.1+0.005*(nstart.p.ros["L"]-38) # make juvenile prawns cost a function of their size with reference at $0.1/juvenile

#post-process to estimate additional parameters
# 0) prawn density rather than raw prawn numbers
  op.ros$P.dens = op.ros$P / 10000
# 1) mean prawn biomass (allometric function)
  op.ros$B = ((par.aqua['a.p']*(op.ros$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op.ros$Bt = op.ros$B*op.ros$P /1000 
# 2.5) marketable prawns at harvest
  eta.ros = predict(eta.lm, newdata = data.frame(dens = nstart.p.ros['P']/area)) # Fraction of harvest that's marketable as function of stocking density
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op.ros$profit = eta.ros*p*(op.ros$Bt/1000)*exp(-delta*(op.ros$t)) - c*(nstart.p["P"]/1000)
# 4) starting total biomass
  start.mass.kg.ros = op.ros$Bt[op.ros$time==0]/1000
# 5) harvest mass in kg (harvest assumed to occur when biomass is maximized)   
  harvest.mass.kg.ros = eta.ros*max(op.ros$Bt)/1000
# 6) average mass of prawns at harvest  
  harvest.size.ros = op.ros$B[op.ros$Bt==max(op.ros$Bt)]
# 7) time of harvest (when biomass is maximized)
  harvest.time.ros = op.ros$time[op.ros$Bt==max(op.ros$Bt)]
# 8) time in months
  op.ros$time.mo = op.ros$time / 30
# 9) species
  op.ros$Species = 'M. rosenbergii'
  
  
#plot showing length trajectory in each species ###########
require(ggplot2)
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
            scale_y_continuous(breaks = c(0, 25, 100, 150,200),
                               labels = c('0', '25', '100', '150', ' 200'),
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
      scale_y_continuous(breaks = seq(0,2, 0.5),
                         labels = c('0', '0.5', '1.0', '1.5', '2.0'),
                         limits = c(0, 2.1))
    
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
    scale_y_continuous(breaks = seq(0,500,100),
                       labels = c('0','100','200','300','400','500'),
                       limits = c(0, 550))
    
    
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
  
  
    



############################


pa.p1 = ggplot(op.vol, aes(x = time)) +
          theme_bw() +
          theme(legend.position = c(0.9, 0.5),        #place legend inside plot
                axis.text = element_text(size = 12),  #increase axis label size
                axis.title = element_text(size = 15), #increase axis title size
                title = element_text(size = 15),      #increase title size
                legend.text = element_text(size = 12),#increase legend text size
                legend.title = element_blank())  +    #suppress legend title
          geom_line(aes(y = P, lty = 'Prawns'), size = 1.25) +
          geom_line(aes(y = Bt/10, lty = 'Biomass'), size = 1.25) +
          labs(x = 'time (days)', y = 'P', lty = "Parameter") +
          scale_y_continuous(sec.axis = sec_axis(~.*.01, name = "Total Prawn Biomass (kg)")) +
          scale_x_continuous(name = 'time (days)', 
                           breaks = sort(c(0, 200, op.vol$time[which(op.vol$Bt == max(op.vol$Bt))], 400, 600)),
                           labels = c('0', '200', expression(italic('T')), '400', '600'),
                           limits = c(0, max(t.p))) +
          geom_vline(xintercept = op.vol$time[which(op.vol$Bt == max(op.vol$Bt))], lty = 3)

pa.p1


#Run model over two years incorporating harvest and restock ###########
  stock = data.frame(var = c('P', 'L'),
                     time = c(harvest.time, harvest.time),
                     value = c(nstart.p['P'], nstart.p['L']),
                     method = rep('rep', 2))
  
op.vol.harvest = as.data.frame(ode(nstart.p,t.p,prawn_biomass,par.aqua, 
                                events = list(data = stock)))
#post-process to estimate additional parameters
# 1) mean prawn biomass (allometric function)
  op.vol.harvest$B = ((par.aqua['a.p']*(op.vol.harvest$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op.vol.harvest$Bt = op.vol.harvest$B*op.vol.harvest$P  
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op.vol.harvest$profit = p*(op.vol.harvest$Bt/1000)*exp(-delta*(op.vol.harvest$t)) - c*(nstart.p["P"]/1000)

#plot
pa.p2 = ggplot(op.vol.harvest, aes(x = time)) +
          theme_bw() +
          geom_line(aes(y = P, lty = 'Prawns'), size = 1.25) +
          geom_line(aes(y = Bt/10, lty = 'Biomass'), size = 1.25) +
          labs(x = 'time (days)', y = 'P',
               title = 'Prawn aquaculture dynamics with harvest at t = 352 days', lty = "Parameter") +
          scale_y_continuous(sec.axis = sec_axis(~.*.01, name = "Total Prawn Biomass (kg)")) +
          theme(legend.position = c(0.9, 0.5),        #place legend inside plot
                axis.text = element_text(size = 12),  #increase axis label size
                axis.title = element_text(size = 15), #increase axis title size
                title = element_text(size = 15),      #increase title size
                legend.text = element_text(size = 12),#increase legend text size
                legend.title = element_blank())       #suppress legend title
pa.p2

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
