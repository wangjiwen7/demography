## Project: Mortality in G7 countries
# G7: United States (US), Canada, France, Germany, Italy, Japan, and the United Kingdom (UK)
## Author: Karen Cheung, Jiwen Wang
## Updated on 2020 Agust 17th

## An R routine for computing M estimates from P-spline-smoothed mortality curves

##To smooth death counts with P-splines over a given ranges of ages, the“Mort1Dsmooth” function of the MortalitySmooth package in R (Camarda 2012) is used.

library("MortalitySmooth")
fit <- "Mort1Dsmooth"(x = age, y = D, offset = log(E))

##where age is a vector with the age values for smoothing (age>=10 in our case), and the vectors D and E denotes the observed death counts and person-years lived in each age interval. 
 #The output of the“Mort1Dsmooth” function is a fit object that contains information about the fitting of the Poisson regression model and about the original data.
 #In order to retrieve our estimation process, and to compute standard errors and hence 95% confidence interval for response and log-rates:
  
summary(fit)
fit1Dfor <- predict(fit, newdata = x,se.fit = TRUE)

se<-fit1Dfor[["se.fit"]]

##The “predict.Mort1Dsmooth” function can be used to obtain predictions of the smoothed mortality curve over very narrow age intervals (i.e. narrower than the original 1-year intervals), 
 #in order to come closer to the concept of force of mortality (or instantaneous death rate). In our case, each 1-year interval is subdivided into 1,000 equal parts, then:
  
delta <- 0.001
age.narrow <- seq(from = min(age), to = max(age), by = delta)
log.mu.hat <- predict(object = fit, newdata = age.narrow)
mu.hat <- exp(log.mu.hat)

##Slope b using age 70 and 90 ,and 80 and 100, is then returned by:

logit.mu.hat <- log(mu.hat/(1-mu.hat))
b1 <- (logit.mu.hat[ age.narrow == 90]-logit.mu.hat[ age.narrow==70])/20
b2 <- (logit.mu.hat[ age.narrow ==100]-logit.mu.hat[ age.narrow==80])/20

##Basic numerical integration methods are used, such as the left Riemann sum, to compute the smoothed survival function corresponding to the P-splines smoothed mortality curve, mu.hat:
  
l.hat <- exp(-cumsum(mu.hat ∗ delta))

##The corresponding smoothed density function, which describes the age-at-death distribution, is then readily obtained:
  
d.hat <- mu.hat ∗ l.hat

##Finally, the estimated value of M is given by the age at which the peak of the heap of deaths occurs in the smoothed density function:
  
M.hat <- age.narrow[which.max(d.hat)]

##Having determined M, its derivatives is then achieved using the data at and after M by:

dM     <- d.hat[age.narrow==M.hat]
d(M+)  <- sum(d.hat)
e(M)   <- sum((age.narrow-M.hat)* d.hat*0.01)/sum(d.hat*0.01)
SD(M+) <- (sum((age.narrow-M.hat)^2* d.hat*0.01)/sum(d.hat*0.01))^0.5
