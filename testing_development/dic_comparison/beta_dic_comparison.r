#dic with beta regression comparison.
#clear environment, load packages.
rm(list=ls())
library(runjags)

#### generate pseudo data. ####
#2 informative covariates, an informative intercept, 1 uninformative covariate.
set.seed(1234)
n.obs <- 200
x1 <- runif(n.obs,20,40)
#x2 <- runif(n.obs, 2, 4)
x3 <- runif(n.obs,0,300)
#x4 <- runif(n.obs,100,110)
intercept <- rep(1, n.obs)
#y1 <- x1*5 + x2*3 + 5 + rnorm(n.obs)
y1 <- x1*5 + 5 + rnorm(n.obs)
y2 <- rep(10, length(y1))
y <- as.matrix(data.frame(y1,y2))
y <- y / rowSums(y)
#x <- as.matrix(data.frame(intercept,x1,x2,x3))
x <- as.matrix(data.frame(intercept,x1,x3))
#specify jags reps throughout.
n.adapt  <- 500
n.burnin <- 1000
n.sample <- 1000

#### Specify a simple beta model. ####
jags.model <- "
model{
  #priors
  for (i in 1:N.preds) {
    m[i] ~ dnorm(0, 1.0E-3)
  }
  tau ~ dgamma(.1,.1)

  #beta part
  for (i in 1:N){
  y[i] ~ dbeta(p[i], q[i])
  p[i] <- mu[i] * tau
  q[i] <- (1 - mu[i]) * tau
  logit(mu[i]) <- inprod(x[i,],m)
  }
} #end mode loop.
"
#### Specify fixed model. No parameter priors because parameters are fixed. ####
fix.model <- "
model{
  #beta part
  for (i in 1:N){
  y[i] ~ dbeta(p[i], q[i])
  p[i] <- mu[i] * tau
  q[i] <- (1 - mu[i]) * tau
  logit(mu[i]) <- inprod(x[i,],m)
  }
} #end model loop
"

#### Fit JAGS model, get dic. ####
jags.data <- list(y=y[,1],x=x,N.preds = ncol(x), N = nrow(y))
jags.out <- run.jags(model = jags.model, data = jags.data, adapt=n.adapt,burnin=n.burnin,sample=n.sample,n.chains=3, monitor = c('m','tau','deviance'))
jags.sum <- summary(jags.out)
jags.dev <- jags.sum[grep('deviance', rownames(jags.sum)),][4]
jags.m   <- jags.sum[grep('m'       , rownames(jags.sum)),][,4]
jags.tau <- jags.sum[grep('tau'     , rownames(jags.sum)),][4]
#get pD with fix.model.
pD.data <- list(y=y[,1],x=x,N.preds = ncol(x), N = length(y[,1]), tau = jags.tau, m = jags.m)
pD.model <- run.jags(model = fix.model, data = pD.data, adapt=100,burnin=100,sample=100,n.chains=3, monitor = c('deviance'))
pD.sum   <- summary(pD.model)
pD       <- pD.sum[3]
dic      <- 2*jags.dev - pD
dic2     <- extract(jags.out, 'dic')
#### leave one predictor out at a time, calculate dic. ####
jags.compare.dic   <- list()
jags.compare.dic2  <- list()
jags.summary.out   <- list()
for(i in 2:ncol(x)){
  new.x <- x[,-i]
  jd <- list(y=y[,1],x=new.x,N.preds = ncol(new.x), N = length(y[,1]))
  jags.out <- run.jags(model = jags.model, data = jd, adapt=n.adapt,burnin=n.burnin,sample=n.sample,n.chains=3, monitor = c('m','tau','deviance'))
  jags.sum <- summary(jags.out)
  jags.dev <- jags.sum[grep('deviance', rownames(jags.sum)),][ 4]
  jags.m   <- jags.sum[grep('m'       , rownames(jags.sum)),][,4]
  jags.tau <- jags.sum[grep('tau'     , rownames(jags.sum)),][ 4]
  #get pD with fix model.
  pD.data <- list(y=y[,1],x=new.x,N = length(y[,1]), tau = jags.tau, m = jags.m)
  pD.model <- run.jags(model = fix.model, data = pD.data, adapt=100,burnin=100,sample=100,n.chains=3, monitor = c('deviance'))
  pD.sum   <- summary(pD.model)
  pD       <- pD.sum[3]
  jags.compare.dic [[i-1]] <- 2*jags.dev - pD
  jags.compare.dic2[[i-1]] <- extract(jags.out, 'dic')
  jags.summary.out [[i-1]] <- jags.sum
}
jags.compare.dic <- unlist(jags.compare.dic)
names(jags.compare.dic) <- colnames(x[,2:ncol(x)])

#### which model improves DIC the most? ####
#two metrics are in line. Both strongly discourage dropping x1 (good).
#both mildly encourage (1-3 points of support) dropping x3 (uninformative predictor).
dic
jags.compare.dic
dic2
jags.compare.dic2
