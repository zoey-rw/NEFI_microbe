#testing dic model selection in JAGS using a linear model.
#This seems to work.
rm(list=ls())
library(runjags)

#### pseudo data generation ####
set.seed(1234)
n <- 200
x1 <- rnorm(n)
x2 <- rnorm(n)
x3 <- rnorm(n)
y <- x1*1.2 + x2*-0.7 + rnorm(n)
summary(lm(y ~ x1 + x2 + x3))
x <- as.matrix(data.frame(x1,x2,x3))

#### Specify JAGS model. ####
jags.model <- "
model {
  #parameter priors
  for(i in 1:n.X){
    m[i] ~ dnorm(0, 1E-3)
  }
  #uncertain prior
    tau <- pow(sigma, -2)
  sigma  ~ dunif(0, 100)
  
  for(i in 1:N){
        y[i] ~ dnorm(y.hat[i], tau)
    y.hat[i] <- inprod(x[i,],m)
  }
}
"
#### Fit JAGS model. ####
jags.data <- list(N=nrow(x), n.X = ncol(x), x=x, y=y)

jags.out <- run.jags(   model = jags.model,
                         data = jags.data,
                        adapt = 100,
                       burnin = 200,
                       sample = 300,
                     n.chains = 3,
                     monitor = c('m'))
dic <- extract(jags.out, 'dic')


#### Fit models without one predictor at a time. ####
#drop each predictor.
jags.compare.data  <- list()
jags.compare.model <- list()
jags.compare.dic   <- list()
for(i in 1:ncol(x)){
  jags.x <- x[,-i]
  jags.compare.data [[i]] <- list(N=nrow(jags.x), n.X = ncol(jags.x), x=jags.x, y=y)
  jags.compare.model[[i]] <- run.jags(model=jags.model, data = jags.compare.data[[i]],adapt=100,burnin=200,sample=300,n.chains=3,monitor=c('m'))
  jags.compare.dic  [[i]] <- extract(jags.compare.model[[i]],'dic')
  
}
#Says to drop 3rd predictor by like 2 points.. Great. Drop and rerun.
dic
jags.compare.dic

#### modify model - inprod doesnt work with just one predictor in the predictor matrix. ####
jags.model <- "
model {
  #parameter priors
  for(i in 1:n.X){
    m[i] ~ dnorm(0, 1E-3)
  }
  #uncertain prior
    tau <- pow(sigma, -2)
  sigma  ~ dunif(0, 100)
  
  for(i in 1:N){
        y[i] ~ dnorm(y.hat[i], tau)
    y.hat[i] <- x[i] * m
  }
}
"
#### try dropping next two predictors which are informative ####
new.x <- x[,-3] #3 was remvoed based on last dic filtering ciriteria.
jags.data <- list(N=nrow(new.x), n.X = ncol(new.x), x=new.x, y=y)
jags.out  <- run.jags(model=jags.model, data = jags.data,adapt=100,burnin=200,sample=300,n.chains=3,monitor=c('m'))
dic <- extract(jags.out, 'dic')
jags.compare.data  <- list()
jags.compare.model <- list()
jags.compare.dic   <- list()
for(i in 1:ncol(new.x)){
  jags.x <- new.x[,-i]
  jags.compare.data [[i]] <- list(N=length(y), n.X = 1, x=jags.x, y=y)
  jags.compare.model[[i]] <- run.jags(model=jags.model, data = jags.compare.data[[i]],adapt=100,burnin=200,sample=300,n.chains=3,monitor=c('m'))
  jags.compare.dic  [[i]] <- extract(jags.compare.model[[i]],'dic')
}

#dic says keep the other predictors for sure.
dic
jags.compare.dic
