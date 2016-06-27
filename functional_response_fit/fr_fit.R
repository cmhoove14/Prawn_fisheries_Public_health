## Fitting alpha and Th functions based on data from Sokolow et al. 2014
## Data extracted from plots using WebPlotDigitizer

a_points = read.csv("functional_response_fit/a_points.csv", header = FALSE)
th_points = read.csv("functional_response_fit/th_points.csv", header = FALSE)

a_model = lm(a_points[,2] ~ 0 + a_points[,1])
summary(a_model)
plot(a_points[,1], a_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Attack rate")
lines(a_points[,1], 0.037192*a_points[,1])

th_model = nls(V2 ~ 1/(c*V1), data = th_points, start = list(c=1))
summary(th_model)
x = c(0:250)
plot(th_points[,1], th_points[,2], xlab = "Prawn-to-snail biomass ratio", ylab = "Handling time")
lines(x, 1/(0.40450*x))

