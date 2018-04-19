#Worldclim2 temperature model.
#clearn environment, load packages.
rm(list=ls())
library(runjags)

#output path for temperature summary object.
out.path <- '/fs/data3/caverill/NEFI_microbial/worldclim2_uncertainty/temperature_JAGS_model.rds'

#load data
temp <- readRDS('/fs/data3/caverill/NEFI_microbial/worldclim2_uncertainty/worldclim2_tavg_uncertainty.rds')

#subset prec for prelim fitting.
temp <- temp[sample(nrow(temp),30000),]



#setup runjags model. Temperature has constant variance.
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

#Setup data object.
temp.data <- list(      y = temp$observed,
                        x = temp$predicted,
                        N = length(temp$observed),
                        N.preds = 2)

#pick starting values based on lm fits to help chain mixing and speed convergence.
start.temp <- list()
temp.freq <- lm(observed ~ predicted, data = temp)
start.temp[[1]] <- list(sigma = 1, m = c(coef(temp.freq)))
start.temp[[2]] <- lapply(start.temp[[1]],'*',1.01)
start.temp[[3]] <- lapply(start.temp[[1]],'*',0.99)

#fit temperature model. This number of samples gives clean looking distributions.
temp.out <- run.jags(temp.model,
                     data = temp.data,
                     adapt = 200,
                     burnin = 1000,
                     sample = 3000,
                     n.chains = 3,
                     monitor = c('m','sigma'),
                     inits = start.temp
)
#summarize and plot output.
temp.summary <- data.frame(summary(temp.out))
plot(temp.out)

#save output.
saveRDS(temp.summary,out.path)