# Data from Ranjeet and Kurup 2006 extracted from table 1 #######

rk06 = data.frame(dens_ha = c(14000, 25000, 40000, 60000),      # Stocking density per hectare (P_nought)    
                  dens = c(14000, 25000, 40000, 60000)/10000,   # Stocking density per m^2
                  surv = c(0.42, 0.38, 0.34, 0.23),             # proportion survival of all stocked prawns
                  marketable = c(0.58, 0.48, 0.43, 0.24),       # prportion of harvested prawns that were marketable
                  mean_B = c(56.8, 51.3, 41.2, 36.7),           # mean prawn mass in grams at harvest  
                  mean_B_perday = c(56.8, 51.3, 41.2, 36.7)/(8*30), # mean daily growth (in grams) 
                  harvest = c(320, 480, 630, 510),              # total prawns harvested in kg
                  market_harvest = c(278, 314, 460, 412),       # marketable prawns harvested in kg
                  profit = c(60166, 87216, 114408, 86241)/48.89)# profit per hectare converted from rupees (reported) to USD (48.89 ruppes/USD reported conversion)
