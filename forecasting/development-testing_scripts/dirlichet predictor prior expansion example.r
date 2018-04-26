rm(list=ls())
#generate 3 columns of species proprotional abundance data
y <- matrix(ncol = 3, nrow = 100)
y[,] <- abs(rnorm(length(y)))
for(i in 1:nrow(y)){
  y[i,] <- y[i,] / sum(y[i,])
}

#generate 2 columns of predictors and an intercept
x <- matrix(ncol = 2, nrow = 100)
x[,] <- rnorm(length(x), mean = 20, sd = 4)
x <- cbind(rep(1,nrow(x)),x)
#x <- data.frame(x)
#colnames(x) <- c('x1','x2','x3')

jags.model = "
model {
    #setup parameter priors for each species * predictor combination.
    for(j in 1:N.spp){
      for(k in 1:N.preds){
        m[k,j] ~ dgamma(1.0E-3, 1.0E-3)
      }
    }

    #go ahead and fit means of species abundances as a linear combination of predictor and parameters.
    for(i in 1:N){
        for(j in 1:N.spp){
             log(a0[i,j]) <- inprod(m[,j] , x[i,])
           }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
    }

} #close model loop.
"
jags.model = "
model {
  #setup parameter priors for each species * predictor combination.
  for(j in 1:N.spp){
      for(k in 1:N.preds){
        m[k,j] ~ dgamma(1.0E-3, 1.0E-3)
      }
  }
  
  #go ahead and fit means of species abundances as a linear combination of predictor and parameters.
  for(i in 1:N){
    for(j in 1:N.spp){
      log(a0[i,j]) <- inprod(m[,j] , x[i,])
    }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }
  
} #close model loop.
"
jags.data <- list(y = as.matrix(y), x = as.matrix(x),
                  N.spp = ncol(y), N.preds = ncol(x), N = nrow(y))


jags.out <- run.jags(jags.model,
                     data=jags.data,
                     adapt = 100,
                     burnin = 200,
                     sample = 400,
                     n.chains=3,
                     monitor=c('m'))
