#clear R environment, load runjags.
rm(list=ls())
library(runjags)

#Simulate some species count data.
y1 <- round(rnorm(100, 100, 5)) 
y2 <- round(y1*-2 + rnorm(100,0,5)) + 400
y3 <- round(rnorm(100,50,5))
spp.y <- data.frame(y1,y2,y3)
spp.y <- as.matrix(spp.y)

#Specify multivariate normal with covariance.
jags.model = "
model {
#setup priors for each species
  for(j in 1:N.spp){
    m0[j] ~ dnorm(0, .001) #intercept prior
  }
  
  #Sigma (covariance matrix) prior
  sigma ~ dwish(r,S)
  
  #implement multivariate normal.
  for(i in 1:N){
    y.hat[i,1:N.spp] <- (m0[1:N.spp])
        y[i,1:N.spp] ~ dmnorm(y.hat[i,1:N.spp], sigma[1:N.spp,1:N.spp])
  }

  #predict species abundances.
#  for(i in 1:N.pred){
#    y.pred[i,1:N.spp] ~ dmnorm(y.hat[i,1:N.spp], sigma[1:N.spp,1:N.spp])
#  }
  

} #close model loop.
"

#Wish prior on covariance - weak covariance.
W = diag(ncol(spp.y))/10000

#get a prediction prior that varies in the abundance of y1
y.pred <- as.matrix(data.frame(round(rnorm(100, 100, 10)), NA, NA))
colnames(y.pred) <- NULL

#Setup JAGS data object.
jags.data <- list(y = spp.y,
                  N = nrow(spp.y),
                  N.spp = ncol(spp.y),
                  r = W,
                  S = ncol(spp.y),
                  N.pred = 100)

#Fit JAGS model.
jags.out <- run.jags(jags.model,
                     data=jags.data,
                     adapt = 100,
                     burnin = 3000,
                     sample = 2000,
                     n.chains=3,
                     monitor=c('m0','sigma'))
j.sum <- summary(jags.out)

#check that intercepts and species mean abundances are the same. They are.
summary(jags.out)[grep('m0',rownames(summary(jags.out))),'Median']
colMeans(spp.y)

#Check that true covariance matrix and estimated covariance matrix are the same.
z <- summary(jags.out)[grep('sigma',rownames(summary(jags.out))),'Median']
i.covm <- matrix(z, ncol = sqrt(length(z)))
solve(i.covm) #jags returns the inverse of the covariance matrix. This inverts it.
cov(spp.y) #true covariance matrix.

