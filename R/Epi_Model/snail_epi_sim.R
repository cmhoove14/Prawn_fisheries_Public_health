sim_epi_mod <- function(pars, nstart = allvh.eqbm.imm){
  
  run_eqbm <- as.data.frame(ode(nstart,
                                seq(1,365*100,50),
                                snail_epi_allvh_imm,
                                pars,
                                atol = 1e-6, rtol = 1e-6, method = "radau")) %>% 
    mutate(W = pars["cvrg"]*Wt + (1-pars["cvrg"])*Wu,
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t)) # density snails of size class 3
  
  return(run_eqbm %>% filter(time == max(time)) %>% select(W:N.t))
}