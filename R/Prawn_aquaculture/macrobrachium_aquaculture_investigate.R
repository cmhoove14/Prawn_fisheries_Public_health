# Data from various aquaculture studies compiled for use in prawn aqauclture model #############
  require(ggplot2)

source("Data/Ranjeet_Kurup_06_data.R")
source("Data/aquaculture_parameters.R")

  rk06$p_kg = rk06$profit / rk06$market_harvest                 # estimate of price: USD per kg marketable harvested

  ggplot(rk06, aes(x = dens, y = marketable)) + 
    theme_bw() +
    labs(x = 'stocking density (P/m^2)', y = 'proportion hravest marketable',
         title = 'Proportion of harvest that is marketable as function of stocking density \n Ranjeet & Kurup 2006 data') +
    xlim(0,7.5) +
    ylim(0,1) +
    geom_point() +
    stat_smooth(method = "lm", col = "darkgreen")
  
  eta.lm = lm(marketable ~ dens, data = rk06)
  
om.gg = ggplot(rk06, aes(x = dens, y = surv)) + 
                theme_bw() +
                labs(x = 'stocking density (P/m^2)', y = '% Survival',
                     title = 'Cohort survival as function of stocking density \n Ranjeet & Kurup 2006 data') +
                xlim(0,7.5) +
                ylim(0,1) +
                geom_point() +
                stat_smooth(method = "lm", col = "red")
  
om.gg

  om.lm = lm(surv ~ dens, data = rk06)
  
gam.gg = ggplot(rk06, aes(x = dens, y = mean_B)) + 
                  theme_bw() +
                  labs(x = 'stocking density (P/m^2)', y = 'Mean size (g)',
                       title = 'Mean prawn weight at harvest as function of stocking density \n Ranjeet & Kurup 2006 data') +
                  xlim(0,7.5) +
                  #ylim(0,1) +
                  geom_point() +
                  stat_smooth(method = "lm", col = "blue")
gam.gg  

  gam.lm = lm(mean_B_perday ~ dens, data = rk06)

p.gg = ggplot(rk06, aes(x = dens, y = profit)) + 
            theme_bw() +
            labs(x = 'stocking density (P/m^2)', y = 'Profit (USD)',
                 title = 'Profit from harvest as function of stocking density \n Ranjeet & Kurup 2006 data') +
            xlim(0,7.5) +
            #ylim(0,1) +
            geom_point()
  p.gg  
  

  
# Lalringsanga 2012 L-W relationship #####
f1 <- function(l, a, b){
  10^a*(l/10)^b  #convert from mm to cm by dividing by 10
}  
  
l_vec <- c(1:300)   

plot(l_vec, sapply(l_vec, f1, a = -2.0570, b = 2.9506), type = 'l', lwd = 2,  #Relationship for juveniles
     xlab = "length(mm)", ylab = "weight (g)", ylim = c(0, 300), xlim = c(0,300))
  lines(l_vec, sapply(l_vec, f1, a = -2.3344, b = 3.2944), lwd = 2, col = 2)  # growout phase
  lines(l_vec, sapply(l_vec, f1, a = -2.2589, b = 3.2667), lwd = 2, col = 3)  # broodstock phase
  lines(l_vec, sapply(l_vec, f1, a = -2.6132, b = 3.5502), lwd = 2, col = 4)  # all males
  lines(l_vec, sapply(l_vec, f1, a = -2.4339, b = 3.3893), lwd = 2, col = 5)  # all pooled

#Also look at data from Kuris et al 1987
  
f2 <- function(l, a, b){
  exp((log(l/a)/b))
}  

lines(l_vec, sapply(l_vec, f2, a = 37, b = 0.298), lwd = 2, col = 6)  # kuris data
