prawn_sim <- function(t.runup, t.intervene, t.post,
                      model.runup = snail_epi_allvh_imm,
                      model.intervene = snail_prawn_model_imm,
                      model.post = snail_epi_allvh_imm,
                      nstart, parameters, stock.events){
  
# Simulate model for runup time ###########   
  runup <- as.data.frame(ode(nstart,
                             t.runup,
                             model.runup,
                             parameters)) %>% 
    mutate(P = 0,
           L = 0,
           W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = par.snails.imm['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3  
  
  nstart.intervene <- runup %>% filter(time == max(t.runup)) %>% select(S1:L) %>% unlist()
    nstart.intervene["P"] <- stocks.ros %>% filter(var == "P" & time == 366) %>% pull(value)
    nstart.intervene["L"] <- stocks.ros %>% filter(var == "L" & time == 366) %>% pull(value)
 
# Simulate model for intervention time ###########        
  intervene <- as.data.frame(ode(nstart.intervene,
                                 t.intervene,
                                 model.intervene,
                                 parameters,
                                 events = list(data = stock.events))) %>% 
    mutate(W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = par.snails.imm['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3
    
  nstart.post <- intervene %>% filter(time == max(t.intervene)) %>% select(S1:Wu) %>% unlist()

# Simulate model post intervention ######
  post <- as.data.frame(ode(nstart.post,
                            t.post,
                            model.post,
                            parameters)) %>% 
    mutate(P = 0,
           L = 0,
           W = cvrg*Wt + (1-cvrg)*Wu,
           prev = pnbinom(2, size = par.snails.imm['phi'], mu = W, lower.tail = FALSE),
           S.t = (S1 + S2 + S3) / area,  # density susceptible snails
           E.t = (E1 + E2 + E3) / area,  # density exposed snails 
           I.t = (I1 + I2 + I3 ) / area, # density infected snails
           N.t = (S.t + E.t + I.t), # total density snails
           t.1 = (S1 + E1 + I1) / area, # density snails of size class 1
           t.2 = (S2 + E2 + I2) / area, # density snails of size class 2
           t.3 = (S3 + E3 + I3) / area) # density snails of size class 3
  
  sim.fin <- rbind(runup[-max(t.runup),], intervene[-which(intervene$time == max(t.intervene)),], post)
  
  return(sim.fin)
}