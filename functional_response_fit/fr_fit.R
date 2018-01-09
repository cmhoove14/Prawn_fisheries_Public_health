## Fitting alpha and Th functions based on data from Sokolow et al. 2014
## Data extracted from plots using WebPlotDigitizer

#Attack rate model
a_points = read.csv("functional_response_fit/a_points.csv", header = FALSE)

x = c(0:2500)

#a_model = nls(V2 ~ -log(3) + b*log(V1), data = a_points, start = list(b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
#summary(a_model)
#lines(x, -log(3) + a_model$m$getPars()*log(x))

plot(a_points[,1], a_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Attack rate", pch = 16,
     xlim = c(0,500))

#logarithmic function
a_model = nls(V2 ~ b*log(V1), data = a_points, start = list(b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model)
  lines(x, a_model$m$getPars()[1]*log(x))
  
#logarithmic function
a_model = nls(V2 ~ c + b*log(V1), data = a_points, start = list(c = -log(3), b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model)
  lines(x, a_model$m$getPars()[1] + a_model$m$getPars()[2]*log(x))

#logarithmic model excluding points that fall in refuge biomass ratio (~6)
a2 = subset(a_points, V1 >= 6)  
  a_model2 = nls(V2 ~ b*log(V1) - c, data = a2, start = list(b=1, c = 3))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model2)
  lines(x, a_model2$m$getPars()[1]*log(x) - a_model2$m$getPars()[2])

#linear model
a_model3 = nls(V2 ~ a+b*V1, data = a_points, start = list(a = 0, b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
  summary(a_model3)
  lines(x, a_model3$m$getPars()[1] + a_model3$m$getPars()[2]*x)
  
AIC(a_model, a_model2, a_model3)  



#hadnling time model  
th_points = read.csv("functional_response_fit/th_points.csv", header = FALSE)

plot(th_points[,1], th_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Handling time", pch = 16)

th_model = nls(V2 ~ 1/(c*V1), data = th_points, start = list(c=1))
summary(th_model)
lines(x, 1/(th_model$m$getPars()*x))

