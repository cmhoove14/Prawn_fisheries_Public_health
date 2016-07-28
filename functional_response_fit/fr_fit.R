## Fitting alpha and Th functions based on data from Sokolow et al. 2014
## Data extracted from plots using WebPlotDigitizer

a_points = read.csv("functional_response_fit/a_points.csv", header = FALSE)
th_points = read.csv("functional_response_fit/th_points.csv", header = FALSE)

x = c(0:250)

a_model = nls(V2 ~ -log(3) + b*log(V1), data = a_points, start = list(b=1))  # Alternate model (linear with forced 0-intercept): lm(a_points[,2] ~ 0 + a_points[,1])
summary(a_model)
plot(a_points[,1], a_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Attack rate")
lines(x, -log(3) + a_model$m$getPars()*log(x))

th_model = nls(V2 ~ 1/(c*V1), data = th_points, start = list(c=1))
summary(th_model)
plot(th_points[,1], th_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Handling time")
lines(x, 1/(th_model$m$getPars()*x))

