#testing calculating DIC from deviance scores vs. dic.samples()
#clear environment, load packages
rm(list=ls())
library(runjags)

#generate a linear relatioship
set.seed(1234)
n <- 100
x <- rnorm(n)
y <- x + rnorm(n)

#make JAGS data object
jags.data <- list(n = n, x =  x, y = y)

#write JAGS model
jags.model = "
model {
  #parameter priors
      m  ~ dnorm(0, 1E-3)
    tau <- pow(sigma, -2)
  sigma  ~ dunif(0, 100)
  
  #model
  for(i in 1:n){
        y[i]  ~ dnorm(y.hat[i], tau)
    y.hat[i] <- x[i]*m
  }
} #end model loop.
"

#model for running with fixed parameter values.
fix.model = "
model {
  #model
  for(i in 1:n){
        y[i]  ~ dnorm(y.hat[i], tau)
    y.hat[i] <- x[i]*m
  }
} #end model loop.
"

jags.out <- run.jags(model = jags.model,
                    data = jags.data,
                    n.chains = 3,
                    adapt = 100,
                    burnin = 100,
                    sample = 1000,
                    monitor = c('m','deviance','tau'))
jags.sum <- summary(jags.out)

#grab mean deiviance and parameter means.
jags.deviance.mu <- jags.sum[2,4]
               m <- jags.sum[1,4]
             tau <- jags.sum[3,4]

#new jags data
jags.data <- list(y=y, x=x, n=n, m=m, tau = tau)

#get pD which is deviance at mean parameter values.
pD.out <- run.jags(model = fix.model,
                   data = jags.data,
                   n.chains = 3,
                   adapt = 100,
                   burnin = 100,
                   sample = 1000,
                   monitor = c('deviance'))
pD.sum <- summary(pD.out)
pD <- pD.sum[1,3]

#Calculate DIC.
dic_A <- 2*jags.deviance.mu - pD

#get deviance using dic.samples
dic_B <- extract(jags.out,'DIC')

#compare calculated DIC scores.
dic_A #deviance calculated by rerunning model with parameters set to their mean values.
dic_B #deviance calculated from dic.samples()