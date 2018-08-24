#propogating worldclim2 uncertainty.
#1. prec model working great now. needed to sample more of the data. uncertainty increses with elevation. DONE.
###temperature- just stop fitting increasing variance with elevation. Fick seems fine with it. DONE.
#2. save these fits with more samples. DONE
#3. save summaries with variable names appended. DONE
#4. simulate variance or sd for an observation, given predicted and elevation. DONE.
#5. convert into a function.
#clear environment, load packages.
rm(list=ls())
library(runjags)

#load data
prec <- readRDS(wc_prec_raw_data.path)
temp <- readRDS(wc_temp_raw_data.path)

#subset prec for prelim fitting.
prec <- prec[sample(nrow(prec),20000),]
temp <- temp[sample(nrow(temp),10000),]

#prec_output.path
#temp_output.path

#plot(observed ~ predicted, data = prec, cex = 0.4, pch = 16, ylim = c(0,4000))

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

#temperature model does not have an uncertainty that scales with elevation.
temp.model = "
model {
  #setup priors for model predictors.
  for(j in 1:N.preds){
    m[j] ~ dnorm(0, .0001) #flat uninformative priors.
  }
  
  #tau is still a function of sigma.
  #sigma is now a funtion of an intercept (as before) and elevation, with a parameter that gets fitted.
  sigma ~ dunif(0,100)
    tau <- pow(sigma, -2)
  
  #fit a linear model
  #should be a way to store predictors as a list with an intercept.
  #then just multiply the m vector by the predictor list for the linear component of the model.
  for(i in 1:N){
    y.hat[i] <- m[1]*1 + m[2]*x[i]
        y[i] ~ dnorm(y.hat[i], tau)
  }
} #close model loop.
"


#Setup data objects.
prec.data <- list(      y = prec$observed,
                        x = prec$predicted,
                     elev = prec$elev,
                        N = length(prec$observed),
                  N.preds = 2)
temp.data <- list(      y = temp$observed,
                        x = temp$predicted,
                        N = length(temp$observed),
                  N.preds = 2)

#pick starting values based on lm fits to help chain mixing and speed convergence.
start.prec <- list()
prec.freq <- lm(observed ~ predicted, data = prec)
start.prec[[1]] <- list(k0 = 47, k1 = 0.003, m = c(coef(prec.freq)))
start.prec[[2]] <- lapply(start.prec[[1]],'*',1.01)
start.prec[[3]] <- lapply(start.prec[[1]],'*',0.99)
start.temp <- list()
temp.freq <- lm(observed ~ predicted, data = temp)
start.temp[[1]] <- list(sigma = 1, m = c(coef(temp.freq)))
start.temp[[2]] <- lapply(start.temp[[1]],'*',1.01)
start.temp[[3]] <- lapply(start.temp[[1]],'*',0.99)

#fit precip model.
prec.out <- run.jags(prec.model,
                     data=prec.data,
                     adapt = 100,
                     burnin = 200,
                     sample = 400,
                     n.chains=3,
                     monitor=c('m','k0','k1'),
                     inits = start.prec)
summary(prec.out)
#plot(prec.out)

#fit temperature model.
temp.out <- run.jags(temp.model,
                       data = temp.data,
                      adapt = 100,
                     burnin = 300,
                     sample = 500,
                   n.chains = 3,
                    monitor = c('m','sigma'),
                      inits = start.temp
                     )
summary(temp.out)
plot(temp.out)

###check variance as a function of elevation###
#generate a single value of precip predicted at 101 different elevations.
d.pred <- data.frame(elevation = seq(0,1780, by = 17.8))
d.pred$predicted <- 1000

#grab parameter means and SD from precip model
parms <- summary(prec.out)[,4:5]

#sample from paramter distributions.
n.samps <- 100
out <- list()
for(i in 1:n.samps){
  #draw parameteres from distributions.
  m1 <- rnorm(1,parms[1,1],parms[1,2])
  m2 <- rnorm(1,parms[2,1],parms[2,2])
  k0 <- rnorm(1,parms[3,1],parms[3,2])
  k1 <- rnorm(1,parms[4,1],parms[4,2])
  
  #esitmate worldclim values based on parameter draws.
     sigma <- d.pred$elevation*k1 + k0
       tau <- sigma^-2
     y.hat <- m1 + m2*d.pred$predicted
         y <- rnorm(length(y.hat),y.hat, sigma)
  out[[i]] <- y
}
#summarize output and plot sd as a function of elevation.
out <- data.frame(do.call('rbind',out))
colnames(out) <- d.pred$elevation
d.pred$out.sd <- apply(out, 2, sd)
plot(d.pred$out.sd ~ d.pred$elevation)
