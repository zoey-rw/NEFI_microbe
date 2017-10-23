#developing a forecasting model in JAGS
#October 19, 2017
#clear enironment, load packages.
rm(list=ls())
library(runjags)
library(boot) #inv.logit function
library(geosphere)

#load prior data and forecast data.
prior.d <- readRDS('forecasting/fungal_prior_data.rds')
  for.d <- readRDS('forecasting/fungal_forecast_data.rds')
  all.d <- rbind(prior.d,for.d)

#generate a temporal distance matrix.
d.time  <- as.matrix(dist(all.d$epoch_date))

#generate a spatial distance matrix.
d.space <- distm(cbind(all.d$longitude, all.d$latitude), fun = distHaversine)

nrow(all.d[all.d$relEM ==1,])

#assign dependent variable
y <- prior.d$relEM
x1 <- prior.d$map30
y.new <- for.d$relEM
x1.new <- for.d$map30

#transform [0,1] to (0,1) a la Cribari-Neto & Zeileis 2010
y <- (y * (length(y) - 1) + 0.5) / length(y)
y.new <- (y.new * (length(y.new) - 1) + 0.5) / length(y.new)


jags.model = "
model{
# priors
a0 ~ dnorm(0, .001)
a1 ~ dnorm(0, .001)
tau ~ dgamma(.1,.1)

# likelihood for mu and tau- This is a beta regression on the continuous values between 0 and 1. mu2 is p, and tau is q. Both are shape parameters. 
for (i in 1:N){
y[i] ~ dbeta(p[i], q[i])
p[i] <- mu[i] * tau
q[i] <- (1 - mu[i]) * tau
logit(mu[i]) <- a0 + a1*x1[i]
}

#forecast loop
for (i in 1:N.new){
y.new[i] ~ dbeta(p.new[i], q.new[i])
p.new[i] <- mu.new[i] * tau
q.new[i] <- (1 - mu.new[i]) * tau
logit(mu.new[i]) <- a0 + a1*x1.new[i]
}



} #end model
"

#list of data
jd <- list(y=y, x1 = x1, x1.new = x1.new, N = length(y), N.new = length(y.new))

#priors
inits <- list()
inits[[1]] <- list(a0 =  0.00)
inits[[2]] <- list(a0 =  0.05)
inits[[3]] <- list(a0 = -0.05)

jags.out <- run.jags(jags.model,
                     data=jd,
                     adapt = 1000,
                     burnin = 1000,
                     sample = 2000,
                     n.chains=3,
                     inits=inits,
                     monitor=c('a0', 'a1'))
summary(jags.out)
plot(jags.out)
