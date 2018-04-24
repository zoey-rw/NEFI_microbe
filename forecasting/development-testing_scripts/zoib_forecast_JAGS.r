#testing a zero one inflated beta regression for relative abundance data from NEON.
rm(list=ls())
library(runjags)

d <- readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/fungal_prior_data.rds')

 y <- d$relEM
x1 <- d$horizon
x2 <- d$doy
x3 <- d$mat30
x4 <- d$map30

# split the data into discrete and continuous components
y.d <- ifelse(y == 1 | y == 0, y, NA)
y.discrete <- ifelse(is.na(y.d), 0, 1)
y.d <- y.d[!is.na(y.d)]
x1.d <- x1[y.discrete == 1]
x2.d <- x2[y.discrete == 1]
n.discrete <- length(y.d)

which.cont <- which(y < 1 & y > 0)
y.c <- ifelse(y < 1 & y > 0, y, NA)
y.c <- y.c[!is.na(y.c)]
n.cont <- length(y.c)
x1.c <- x1[which.cont]
x2.c <- x2[which.cont]
n <- length(y)

jags.model = "
model{
# priors
a0 ~ dnorm(0, .001)
a1 ~ dnorm(0, .001)
b0 ~ dnorm(0, .001)
b1 ~ dnorm(0, .001)
c0 ~ dnorm(0, .001)
c1 ~ dnorm(0, .001)
t0 ~ dnorm(0, .01)
tau <- exp(t0)

# likelihood for alpha - whether the pdf is dominated by discrete or continuous component.
for (i in 1:n){
logit(alpha[i]) <- a0 + a1 * x2[i] 
y.discrete[i] ~ dbern(alpha[i])
}

# likelihood for gamma - probability a discrete value is zero or one.
for (i in 1:n.discrete){
y.d[i] ~ dbern(mu1[i])
logit(mu1[i]) <- b0 + b1 * x2.d[i]
}

# likelihood for continuous component - predicted value on interval (0,1)
for (i in 1:n.cont){
y.c[i] ~ dbeta(p[i], q[i])
p[i] <- mu2[i] * tau
q[i] <- (1 - mu2[i]) * tau
logit(mu2[i]) <- c0 + c1 * x2.c[i]
}

for ( i in 1:n){
#get predictions sending all observations of x through both models.
logit(mu1.all[i]) <- b0 + b1 * x2[i]
logit(mu2.all[i]) <- c0 + c1 * x2[i]
p.all[i] <- mu2.all[i] * tau
q.all[i] <- (1 - mu2.all[i]) * tau

#get a predicted y value - use alpha to blend bernoulli and beta models
bern.weight[i] ~ dbern(alpha[i])
beta.weight[i] ~ dbern(1 - alpha[i])
  bern.comp[i] ~ dbern(mu1.all[i])
  beta.comp[i] ~ dbeta(p.all[i], q.all[i])
     pred.y[i] <- bern.weight[i] * bern.comp[i] + beta.weight[i] * beta.comp[i]
}
}  
"

jd <- list(x2=x2, y.d=y.d, y.c=y.c, y.discrete = y.discrete,
           n.discrete=n.discrete, n.cont = n.cont,
           x2.d = x2.d, x2.c=x2.c, n=n)

jags.out <- run.jags(jags.model,
                     data=jd,
                     adapt = 1000,
                     burnin = 1000,
                     sample = 2000,
                     n.chains=3,
                     #inits=inits,
                     monitor=c('a0', 'a1','b0','b1','c0','c1','tau','pred.y'))

output <- summary(jags.out)

#grab predicted y values
preds <-    data.frame(output[grep('pred.y', rownames(output)),])
pred.y <- preds$Mean

par(mfrow = c(1,3))
plot(y ~ pred.y)
plot(y ~ x2)
plot(pred.y ~ x2)
