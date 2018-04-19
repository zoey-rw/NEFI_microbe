#Worldclim2 precipitation model.
#clearn environment, load packages.
rm(list=ls())
library(runjags)

#output path for precipitation summary object.
out.path <- '/fs/data3/caverill/NEFI_microbial/worldclim2_uncertainty/precipitation_JAGS_model.rds'

#load data
prec <- readRDS('/fs/data3/caverill/NEFI_microbial/worldclim2_uncertainty/worldclim2_prec_uncertainty.rds')

#subset prec for prelim fitting.
prec <- prec[sample(nrow(prec),30000),]



#setup runjags model. precipitation has constant variance.
#setup runjags model. Same model for temperature and precipitation.
prec.model = "
model {
  #setup priors for model predictors.
  for(j in 1:N.preds){
    m[j] ~ dnorm(0, .0001) #flat uninformative priors.
  }
  
  #elevation affects uncertainty.
  #tau is still a function of sigma.
  #sigma is now a funtion of an intercept (as before) and elevation, with a parameter that gets fitted.
  k0 ~ dunif(0,100)
  k1 ~ dnorm(0, .0001)
  
  #fit a linear model
  #should be a way to store predictors as a list with an intercept.
  #then just multiply the m vector by the predictor list for the linear component of the model.
  for(i in 1:N){
    y.hat[i] <- m[1]*1 + m[2]*x[i]
    sigma[i] <- k0 + k1*elev[i]
      tau[i] <- pow(sigma[i], -2)
        y[i] ~ dnorm(y.hat[i], tau[i])
  }
} #close model loop.
"

#Setup data objects.
prec.data <- list(      y = prec$observed,
                        x = prec$predicted,
                        elev = prec$elev,
                        N = length(prec$observed),
                        N.preds = 2)

#pick starting values based on lm fits to help chain mixing and speed convergence.
start.prec <- list()
prec.freq <- lm(observed ~ predicted, data = prec)
start.prec[[1]] <- list(k0 = 47, k1 = 0.003, m = c(coef(prec.freq)))
start.prec[[2]] <- lapply(start.prec[[1]],'*',1.01)
start.prec[[3]] <- lapply(start.prec[[1]],'*',0.99)

#fit precipitation model. This number of samples gives clean looking distributions.
prec.out <- run.jags(prec.model,
                     data = prec.data,
                     adapt = 200,
                     burnin = 1000,
                     sample = 3000,
                     n.chains = 3,
                     monitor = c('m','k0','k1'),
                     inits = start.prec
)
#summarize and plot output.
prec.summary <- data.frame(summary(prec.out))
plot(prec.out)

#save output.
saveRDS(prec.summary,out.path)