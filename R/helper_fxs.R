#Helper functions for dynamic model

fx<-function(x, 
             mean.worm, 
             clump){
  alpha<-(mean.worm)/(clump+mean.worm)
  (1-cos(x)) / ( (1+(alpha*cos(x)))^(1+clump) )
}

phi_Wk<-function(W, 
                 phi){
  alpha<-W/(W+phi)
  1-( (1-alpha)^(phi+1) * (integrate(fx, 0, 2*pi, W, phi, stop.on.error = F)$value) /(2*pi)  )
}

get_prev <- function(clump, 
                     burden){
  pnbinom(1, size = clump, mu = burden, lower.tail = FALSE)
}

#Function estimate schisto DALYs
est_dalys <- function(burden,     # mean worm burden of NB
                      clump,      # clumping parameter of NB
                      weights_lo, # disability weight for low infection intensity 
                      weights_hi, # disability weight for high infection intensity >=50 eggs/10mL
                      epmL,       # estimate of eggs/10mL per mated female worm
                      pop,        # human population
                      daily = TRUE){  # estimating DALYs on a daily basis?
  eggs_mL <- rnbinom(pop, size = clump, mu = burden) * phi_Wk(burden, clump) * 0.5 * epmL # estimate of egg burden converted from worm burden
  
  n_lo <- length(eggs_mL[eggs_mL > 0 & eggs_mL < 50])
  n_hi <- length(eggs_mL[eggs_mL >= 50])
  
  if(daily){
    
    dalys_lo <- n_lo * weights_lo/365
    dalys_hi <- n_hi * weights_hi/365
    
  } else {
    
    dalys_lo <- n_lo * weights_lo
    dalys_hi <- n_hi * weights_hi

  }
  
  return(dalys_hi + dalys_lo)
}

#Function to return summary of a vector as median (IQR)
get_sum <- function(vec){
  paste0(round(median(vec, na.rm = T), 2), 
         " (", round(quantile(vec, 0.25, na.rm = T), 2), 
         " - ", round(quantile(vec, 0.75, na.rm = T), 2), ")"  )
} 
