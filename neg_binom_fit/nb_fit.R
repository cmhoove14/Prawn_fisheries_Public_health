## Fitting clumping parameter k for negative binomial distribution,
## based on burden-of-infection data from EPLS

require(MASS)

infection = read.csv("neg_binom_fit/egg_burden.csv")

dist1 = fitdistr(na.omit(infection$Eggs1), "negative binomial")
dist2 = fitdistr(na.omit(infection$Eggs2), "negative binomial")


