## Fitting alpha and Th functions based on data from Sokolow et al. 2014

#Attack rate model
a_points = read.csv("Data/a_points.csv", header = FALSE)

x = c(0:500)

plot(a_points[,1], a_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Attack rate", pch = 16,
     xlim = c(0,500))

#logarithmic function
a_model0 = nls(V2 ~ b*log(V1), data = a_points, start = list(b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model0)
  lines(x, a_model0$m$getPars()[1]*log(x))
  
#logarithmic function
a_model = nls(V2 ~ c + b*log(V1), data = a_points, start = list(c = -log(3), b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model)
  lines(x, a_model$m$getPars()[1] + a_model$m$getPars()[2]*log(x), col = 2)

#logarithmic model excluding points that fall in refuge biomass ratio (~6)
a2 = subset(a_points, V1 >= 3)  
  a_model2 = nls(V2 ~ b*log(V1) - c, data = a2, start = list(b=1, c = 3))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model2)
  lines(x, a_model2$m$getPars()[1]*log(x) - a_model2$m$getPars()[2], col = 3)

#linear model
a_model3 = nls(V2 ~ a+b*V1, data = a_points, start = list(a = 0, b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model3)
  lines(x, a_model3$m$getPars()[1] + a_model3$m$getPars()[2]*x, col = 4)
  
AIC(a_model0, a_model, a_model2, a_model3)  

#Model 2 provides the best fit

#handling time model  
th_points = read.csv("Data/th_points.csv", header = FALSE)

plot(th_points[,1], th_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Handling time", pch = 16)

th_model = nls(V2 ~ 1/(c*V1), data = th_points, start = list(c=1))
summary(th_model)
lines(x, 1/(th_model$m$getPars()*x))

#Test per capita consumption rate given prawn-snail biomass ratio #########
n_dens = c(0:50)

a_rate <- function(bmr, n, ex, pen){
  a = predict(a_model0, newdata = data.frame(V1 = bmr))
  Th = predict(th_model, newdata = data.frame(V1 = bmr))
  
  (a*n^ex)/(1+a*Th*n^ex)/pen
}

plot(n_dens, sapply(n_dens, a_rate, bmr = 100, ex = 1, pen = 10), type = 'l', ylim = c(0,30),
     xlab = "Snail density", ylab = "Snails consumed per prawn per day", 
     main = "Consumption rate per prawn, Holling's II")
  lines(n_dens, sapply(n_dens, a_rate, bmr = 50, ex = 1, pen = 10), col = 2) 
  lines(n_dens, sapply(n_dens, a_rate, bmr = 200, ex = 1, pen = 10), col = 4) 

  
plot(n_dens, sapply(n_dens, a_rate, bmr = 100, ex = 1.5, pen = 10), type = 'l', ylim = c(0,30),
     xlab = "Snail density", ylab = "Snails consumed per prawn per day", 
     main = "Consumption rate per prawn, Holling's III")
  lines(n_dens, sapply(n_dens, a_rate, bmr = 50, ex = 1.5, pen = 10), col = 2) 
  lines(n_dens, sapply(n_dens, a_rate, bmr = 200, ex = 1.5, pen = 10), col = 4) 
   
  
#Test with multiple size classes #########
a_rate2 <- function(n1, n2, n3, bm1, bm2, bm3, ex, pen){
  a1 = predict(a_model0, newdata = data.frame(V1 = bm1))
  Th1 = predict(th_model, newdata = data.frame(V1 = bm1))
  
  c1 <- (a1*n1^ex)/(1+a1*Th1*(n1+n2+n3)^ex)/pen
  
  a2 = predict(a_model0, newdata = data.frame(V1 = bm2))
  Th2 = predict(th_model, newdata = data.frame(V1 = bm2))
  
  c2 <- (a2*n2^ex)/(1+a2*Th2*(n1+n2+n3)^ex)/pen

  
  a3 = predict(a_model0, newdata = data.frame(V1 = bm3))
  Th3 = predict(th_model, newdata = data.frame(V1 = bm3))
  
  c3 <- (a3*n3^ex)/(1+a3*Th3*(n1+n2+n3)^ex)/pen

  return(c(c1, c2, c3))
}

sum(a_rate2(30,10,5,100,100,100,1,10))
a_rate(100, 45,1,10)

# Consumption rates hold if biomass ratio is consistent across different snail classes

#But how should it perform with different biomass ratios and densities
a_rate2(30,10,5,200,100,50,1,10)
sum(a_rate2(30,10,5,200,100,50,1,10))

# Distribute predation between classes logically 
a_rate3 <- function(n1, n2, n3, bm1, bm2, bm3, ex, pen){
  a1 = predict(a_model0, newdata = data.frame(V1 = bm1))
  Th1 = predict(th_model, newdata = data.frame(V1 = bm1))
  a2 = predict(a_model0, newdata = data.frame(V1 = bm2))
  Th2 = predict(th_model, newdata = data.frame(V1 = bm2))
  a3 = predict(a_model0, newdata = data.frame(V1 = bm3))
  Th3 = predict(th_model, newdata = data.frame(V1 = bm3))
  
  c1 <- (a1*n1^ex)/(1+sum(c(a1*Th1*n1,
                            a2*Th2*n2,
                            a3*Th3*n3))^ex)/pen
  
  
  c2 <- (a2*n2^ex)/(1+sum(c(a1*Th1*n1,
                            a2*Th2*n2,
                            a3*Th3*n3))^ex)/pen

  
  
  c3 <- (a3*n3^ex)/(1+sum(c(a1*Th1*n1,
                            a2*Th2*n2,
                            a3*Th3*n3))^ex)/pen

  return(c(c1, c2, c3))
}

a_rate3(30,10,5,200,100,50,1,10)
sum(a_rate3(30,10,5,200,100,50,1,10))
