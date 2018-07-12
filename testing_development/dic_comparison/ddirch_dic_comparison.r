#ddirch covariate selection based on DIC scores.
#clear environment, load packages.
rm(list=ls())
library(runjags)

#### generate pseudo data. ####
#2 informative covariates, an informative intercept, 1 uninformative covariate.
set.seed(1234)
n.obs <- 200
x1 <- runif(n.obs,20,40)
x2 <- runif(n.obs, 2, 4)
x3 <- runif(n.obs,100,110)
#x4 <- runif(n.obs,100,110)
intercept <- rep(1, n.obs)
y1 <- x1*0.5 + x2*3 + 5 + rnorm(n.obs)
y2 <- rep(10, length(y1))
y <- as.matrix(data.frame(y1,y2))
y <- y / rowSums(y)
x <- as.matrix(data.frame(intercept,x1,x2,x3))

#### Define simples ddirch JAGS model. ####
jags.model = "
model {
  #parameter priors for each species.
  alpha ~ dnorm(0, 1.0E-3) 
  for(i in 1:N.preds){
    x.mm[i,1] <- 0
    for (j in 2:N.spp) {
      x.mm[i,j] ~ dnorm(0, 1.0E-3)
    }
  }

  #mean center all predictors (except intercept).
  for(i in 1:N){
    x.center[i,1] <- 1
    for(j in 2:N.preds){
      x.center[i,j] <- x[i,j] - mean(x[,j])
    }
  }
  
  #save mean values for back transforming intercept values.
  for(j in 1:N.preds){x.center.save[j] <- mean(x[,j])}
  
  #fit species abundances as a linear combination of predictors and parameters.
  for(i in 1:N){
    for(j in 1:N.spp){
      log(a0[i,j]) <- alpha + inprod(x.mm[,j], x.center[i,])
    }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }
  
  #map to original parameterization, assuming first column of predictors is intercept.
  for (j in 1:N.spp) {
    x.m[1,j] <- alpha + x.mm[1,j] - inprod(x.mm[2:N.preds,j], x.center.save[2:N.preds])
    for (i in 2:N.preds){
      x.m[i,j] <- x.mm[i,j]
    }
  }
} #close model loop.
"
#### define model for fixed parameters ####
fix.model = "
model{
  #mean center all predictors (except intercept).
  for(i in 1:N){
    x.center[i,1] <- 1
    for(j in 2:N.preds){
      x.center[i,j] <- x[i,j] - mean(x[,j])
    }
  }
  
  #fit species abundances as a linear combination of predictors and parameters.
  for(i in 1:N){
    for(j in 1:N.spp){
      log(a0[i,j]) <- alpha + inprod(x.mm[,j], x.center[i,])
    }
    y[i,1:N.spp] ~ ddirch(a0[i,1:N.spp]) 
  }
} #close model loop.
"

#### fit JAGS model, get deviance and pD ####
jags.data <- list(y=y,x=x,N.preds = ncol(x), N.spp=ncol(y), N = nrow(y))
jags.out <- run.jags(model = jags.model, data = jags.data, adapt=100,burnin=300,sample=300,n.chains=3, monitor = c('x.mm','alpha','deviance'))
jags.sum <- summary(jags.out)
jags.dev <- jags.sum[grep('deviance', rownames(jags.sum)),][4]
jags.mm  <- jags.sum[grep('x.mm'    , rownames(jags.sum)),][,4]
jags.alp <- jags.sum[grep('alpha'   , rownames(jags.sum)),][4]
#convert jags.mm to matrix of mean values
dim(jags.mm) <- c(ncol(x),ncol(y))
#get pD with fix model.
pD.data <- list(y=y,x=x,N.preds = ncol(x), N = nrow(y), N.spp = ncol(y), alpha = jags.alp, x.mm = jags.mm)
pD.model <- run.jags(model = fix.model, data = pD.data, adapt=100,burnin=300,sample=300,n.chains=3, monitor = c('deviance'))
pD.sum   <- summary(pD.model)
pD       <- pD.sum[3]
dic      <- 2*jags.dev - pD

#### leave one predictor out at a time, calculate dic. ####
jags.compare.dic   <- list()
jags.compare.sum   <- list()
for(i in 2:ncol(x)){
  new.x <- x[,-i]
  jd <- list(y=y,x=new.x,N.preds = ncol(new.x), N.spp=ncol(y), N = nrow(y))
  jags.out <- run.jags(model = jags.model, data = jd, adapt=100,burnin=300,sample=300,n.chains=3, monitor = c('x.mm','alpha','deviance'))
  jags.sum <- summary(jags.out)
  jags.dev <- jags.sum[grep('deviance', rownames(jags.sum)),][ 4]
  jags.mm  <- jags.sum[grep('x.mm'    , rownames(jags.sum)),][,4]
  jags.alp <- jags.sum[grep('alpha'   , rownames(jags.sum)),][ 4]
  #convert jags.mm to matrix of mean values
  dim(jags.mm) <- c(ncol(new.x),ncol(y))
  #get pD with fix model.
  pD.data <- list(y=y,x=new.x,N.preds = ncol(new.x), N = nrow(y), N.spp = ncol(y), alpha = jags.alp, x.mm = jags.mm)
  pD.model <- run.jags(model = fix.model, data = pD.data, adapt=100,burnin=300,sample=300,n.chains=3, monitor = c('deviance'))
  pD.sum   <- summary(pD.model)
  pD       <- pD.sum[3]
  jags.compare.dic[[i-1]] <- 2*jags.dev - pD
  jags.compare.sum[[i-1]] <- jags.sum
}
jags.compare.dic <- unlist(jags.compare.dic)
names(jags.compare.dic) <- colnames(x[,2:ncol(x)])

#### which model improves DIC the most? ####
#dic says don't you dare remove x1 or x2. Supports removing x3 by 2-3 points.
dic
jags.compare.dic

#### Testing with function. ####
source('NEFI_functions/covariate_selection_JAGS.r')
test <- covariate_selection_JAGS(y=y,x=x)


#### plotting y vs preditor relationships. ####
#just to prove- plot y as a function x1-x3
par(mfrow = c(ncol(x)-1,ncol(y)))
for(i in 2:4){
  plot(y[,1] ~ x[,i], main = paste0('y1 ~ x',i))
  abline(lm(y[,1] ~ x[,i]), lwd =2 , col = 'green')
  plot(y[,2] ~ x[,i], main = paste0('y2 ~ x',i))
  abline(lm(y[,2] ~ x[,i]), lwd =2 , col = 'green')
}