#Generate plot of aquaculture cycle over two years to display in model summary figure
source('Prawn_aquaculture/prawn_aquaculture_mod.R')

#run model through 2 years ######
op2 = as.data.frame(ode(nstart.p,t.p,prawn_biomass,par.aqua))

#post-process to estimate additional parameters
# 1) mean prawn biomass (allometric function)
  op2$B = ((par.aqua['a.p']*(op2$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op2$Bt = op2$B*op2$P  
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op2$profit = p*(op2$Bt/1000)*exp(-delta*(op2$t)) - c*(nstart.p["P"]/1000)
# 4) starting total biomass
  start.mass.kg = op2$Bt[op2$time==0]/1000
# 5) harvest mass in kg (harvest assumed to occur when biomass is maximized)   
  harvest.mass.kg = max(op2$Bt)/1000
# 6) average mass of prawns at harvest  
  harvest.size = op2$B[op2$Bt==max(op2$Bt)]
# 7) time of harvest (when biomass is maximized)
  harvest.time = op2$time[op2$Bt==max(op2$Bt)]

#plot to show how total biomass (in kg) and number prawns (in 1000s) change over time
require(ggplot2)
  
pa.p1 = ggplot(op2, aes(x = time)) +
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
                           breaks = sort(c(0, 200, op2$time[which(op2$Bt == max(op2$Bt))], 400, 600)),
                           labels = c('0', '200', expression(italic('T')), '400', '600'),
                           limits = c(0, max(t.p))) +
          geom_vline(xintercept = op2$time[which(op2$Bt == max(op2$Bt))], lty = 3)

pa.p1
#Run model over two years incorporating harvest and restock ###########
  stock = data.frame(var = c('P', 'L'),
                     time = c(harvest.time, harvest.time),
                     value = c(nstart.p['P'], nstart.p['L']),
                     method = rep('rep', 2))
  
op2.harvest = as.data.frame(ode(nstart.p,t.p,prawn_biomass,par.aqua, 
                                events = list(data = stock)))
#post-process to estimate additional parameters
# 1) mean prawn biomass (allometric function)
  op2.harvest$B = ((par.aqua['a.p']*(op2.harvest$L/10)^par.aqua['b.p'])/10)
# 2) total prawn biomass (mean biomass * number of prawns)  
  op2.harvest$Bt = op2.harvest$B*op2.harvest$P  
# 3) profit (in terms of revenue (discounted by time since stocking) minus stocking costs )  
  op2.harvest$profit = p*(op2.harvest$Bt/1000)*exp(-delta*(op2.harvest$t)) - c*(nstart.p["P"]/1000)

#plot
pa.p2 = ggplot(op2.harvest, aes(x = time)) +
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