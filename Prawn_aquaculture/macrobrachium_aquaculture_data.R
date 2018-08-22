# Data from various aquaculture studies compiled for use in prawn aqauclture model #############
  require(ggplot2)

# Ranjeet and Kurup 2006 #######

rk06 = data.frame(dens_ha = c(14000, 25000, 40000, 60000),      # Stocking density per hectare (P_nought)    
                  dens = c(14000, 25000, 40000, 60000)/10000,   # Stocking density per m^2
                  surv = c(0.42, 0.38, 0.34, 0.23),             # proportion survival of all stocked prawns
                  marketable = c(0.58, 0.48, 0.43, 0.24),       # prportion of harvested prawns that were marketable
                  mean_B = c(56.8, 51.3, 41.2, 36.7),           # mean prawn mass in grams at harvest  
                  mean_B_perday = c(56.8, 51.3, 41.2, 36.7)/(8*30), # mean daily growth (in grams) 
                  harvest = c(320, 480, 630, 510),              # total prawns harvested in kg
                  market_harvest = c(278, 314, 460, 412),       # marketable prawns harvested in kg
                  profit = c(60166, 87216, 114408, 86241)/48.89)# profit per hectare converted from rupees (reported) to USD (48.89 ruppes/USD reported conversion)

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
  

  
# Lalringsanga 2012 and Kuris 1987 L-W relationship #####
f1 <- function(l, a, b){
  a/10*(l/10)^b  #convert from mm to cm by dividing by 10
}  
  
l_vec <- c(1:300)   

plot(l_vec, sapply(l_vec, f1, a = 0.127836, b = 2.9506), type = 'l', lwd = 2,  #Relationship for juveniles
     xlab = "length(mm)", ylab = "weight (g)", ylim = c(0, 300), xlim = c(0,300))
  lines(l_vec, sapply(l_vec, f1, a = 0.096868, b = 3.2944), lwd = 2, col = 2)  # growout phase
  lines(l_vec, sapply(l_vec, f1, a = 0.101723, b = 3.2667), lwd = 2, col = 3)  # broodstock phase
  lines(l_vec, sapply(l_vec, f1, a = exp(-2.6132), b = 3.5502), lwd = 2, col = 4)  # all males
  lines(l_vec, sapply(l_vec, f1, a = exp(-2.4339), b = 3.3893), lwd = 2, col = 5)  # all pooled

f_mm <- function(l, a, b){
  a*l^b
} 

  lines(l_vec, sapply(l_vec, f_mm, a = 1.21e-6, b = 3.43), lwd = 2, col = 7)

#Also look at data from Kuris et al 1987
  
f2 <- function(l, a, b){
  exp((log(l/a)/b))
}  

lines(l_vec, sapply(l_vec, f2, a = 37, b = 0.298), lwd = 2, col = 6)  # kuris data

# This produces a much shallower (and more reasonable) curve. 
# To get the length-weight rather than weight-length relationship, we'll refit the parameters

plot(log(l_vec), log(sapply(l_vec, f2, a = 37, b = 0.298)),
     xlab = 'log(length)', ylab = 'log(weight)')

kuris_lm <- lm(log(sapply(l_vec, f2, a = 37, b = 0.298)) ~ log(l_vec))
  abline(coef(kuris_lm), lty = 2, col = 2)

kuris_a <- as.numeric(exp(coef(kuris_lm)[1]))
kuris_b <- as.numeric(coef(kuris_lm)[2])

plot(l_vec, sapply(l_vec*10, f1, a = kuris_a*10, b = kuris_b), type = 'l', lwd = 2,  
     xlab = "length(mm)", ylab = "weight (g)", ylim = c(0, 300), xlim = c(0,300))
