#Developing multivariate normal with covariance. This recovers species abundances and covariance matrix.
#Eventually extend to dirlichet case.
#clear ennvironment, load packages.
rm(list=ls())
library(runjags)

#try with simulated data
y1 <- round(rnorm(100, 100, 5)) 
y2 <- round(y1*-2 + rnorm(100,0,5)) + 400
y3 <- round(rnorm(100,50,5))
y4 <- round(rnorm(100,50,5))
spp.y <- data.frame(y1,y2,y3,y4)
spp.y <- as.matrix(spp.y)
seqs <- rowSums(spp.y)

#introduce a sequence depth effect.
#Get a integer value to multiply species abundances by.
#effort <- round(runif(100,5,50))
#spp.y2 <- spp.y * effort
#eqs <- rowSums(spp.y2)
#seqs <- rep(1,nrow(spp.y))


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
  y.hat[i,1:N.spp] <- (m0[1:N.spp]) * seqs[i]
      y[i,1:N.spp] ~ dmnorm(y.hat[i,1:N.spp], sigma[1:N.spp,1:N.spp])
}

} #close model loop.
"

#Wish prior on covariance - weak covariance.
W = diag(ncol(spp.y))/10000

jags.data <- list(y = spp.y,
                  N = nrow(spp.y),
                  N.spp = ncol(spp.y),
                  r = W,
                  S = ncol(spp.y),
                  seqs = seqs)

jags.out <- run.jags(jags.model,
                     data=jags.data,
                     adapt = 100,
                     burnin = 3000,
                     sample = 2000,
                     n.chains=3,
                     monitor=c('m0','sigma'))

#true relative abundances across samples.
colMeans(spp.y) / sum(colMeans(spp.y)) 
#estimated relative abundances
summary(jags.out)[grep('m0',rownames(summary(jags.out))),'Median']


#Get estimated covariance matrix
z <- summary(jags.out)[grep('sigma',rownames(summary(jags.out))),'Median']
z.lo.95 <- summary(jags.out)[grep('sigma',rownames(summary(jags.out))),'Lower95']
z.hi.95 <- summary(jags.out)[grep('sigma',rownames(summary(jags.out))),'Upper95']
i.covm <- matrix(z, ncol = sqrt(length(z)))
i.covm.lo <- matrix(z.lo.95, ncol = sqrt(length(z.lo.95)))
i.covm.hi <- matrix(z.hi.95, ncol = sqrt(length(z.hi.95)))
solve(i.covm) #jags returns the inverse of the covariance matrix. This inverts it.

#What is the covariance of real data?
#spp.y.cov <- spp.y2 / seqs
#cov(spp.y.cov)
cov(spp.y)
cov(spp.y2)

cov2cor(cov(spp.y2))
(cov2cor(cov(spp.y)))^2

cov2cor(solve(i.covm))
cov2cor(solve(i.covm.lo))
cov2cor(solve(i.covm.hi))
