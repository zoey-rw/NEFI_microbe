#Developing dirlichet distribution for EM-SAP abundance.
#clear ennvironment, load packages.
rm(list=ls())
library(runjags)

#load data
d <- readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_all_prior_data.rds')

#specify your y dependent variable matrix of species abundances
#spp.y <- data.frame(d$Ectomycorrhizal,d$Saprotroph)
#get a 3rd species abundance column, abundance of all other fungi that are not ecto or sap.
#spp.y$S3other <- 1 - (rowSums(spp.y))
#names(spp.y) <- c('y1','y2','y3')
#spp.y <- as.matrix(spp.y)
#x1 <- d$pH

#Simulate some species count data.
y1 <- round(rnorm(100, 100, 5)) 
y2 <- round(y1*-2 + rnorm(100,0,5)) + 400
y3 <- round(rnorm(100,50,5))
spp.y <- data.frame(y1,y2,y3)
spp.y <- as.matrix(spp.y)
spp.y <- spp.y/rowSums(spp.y)


jags.model = "
model {
    #setup priors for each species
    for(j in 1:N.spp){
      m0[j] ~ dgamma(1.0E-3, 1.0E-3) #intercept prior
    }

    #Sigma (covariance matrix) prior
    sigma ~ dwish(r,S)

    #implement dirlichet
    for(i in 1:N){
    for(j in 1:N.spp){
#         log(a0[i,j]) <- m0[j] 
       log(a0[i,j]) <- mu[i,j]
     }

     mu[i,1:N.spp] ~ dmnorm(m0,sigma)

     
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
    #y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp], sigma[1:N.spp,1:N.spp]) 
    }

} #close model loop.
"
#Wish prior on covariance - weak covariance.
W = diag(ncol(spp.y))/10000

jags.data <- list(y = spp.y,
                  N = nrow(spp.y),
                  N.spp = ncol(spp.y),
                  r = W,
                  S = ncol(spp.y))

jags.out <- run.jags(jags.model,
                     data=jags.data,
                     adapt = 100,
                     burnin = 200,
                     sample = 400,
                     n.chains=3,
                     monitor=c('m0','mu','sigma'))
summary(jags.out)
mod.sum <- summary(jags.out)

#Can we capture the means? yes.
colMeans(spp.y)
a.m <- exp(summary(jags.out)[,"Median"])
a.m / sum(a.m)

z <- mod.sum[grep('sigma',rownames(mod.sum)),'Median']
i.covm <- matrix(z, ncol = sqrt(length(z)))
solve(i.covm)

cov(spp.y)

cor(spp.y)


jags.mat <- as.matrix(coda::as.mcmc.list(jags.out))
nsim <- 1000
p.sim <- matrix(NA,nsim,ncol(spp.y))
for(i in 1:nsim){
  a0 <- exp(jags.mat[sample.int(nrow(jags.mat),1),])
  p.sim[i,] <- gtools::rdirichlet(1,a0)
}
cor(p.sim)
cor(spp.y)
