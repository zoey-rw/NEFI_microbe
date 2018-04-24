#developing a forecasting model in JAGS - using a beta distribution, will expand to zoib.
#October 19, 2017
#I've successfully implemented a hierarchical structure and random effects.
#It runs but chains do not converge. I've tried picking starting values based on lme fit with no luck.
#clear enironment, load packages.
rm(list=ls())
library(runjags)
library(data.table)
library(parallel)

#specify output paths
o.dir <- '/fs/data3/caverill/NEFI_microbial/prior_data/'
jags.out.path <- paste0(o.dir,'test_jags_out.rds')
data.out.path <- paste0(o.dir,'test_data_out.rds')

#load Talbot 2014 PNAS data
core.d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/tal_core_data.rds'))
site.d <- data.table(readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/tal_site_data.rds'))

#subset to complete cases
core.d <- core.d[,soilTemp := NULL]
core.d <- core.d[complete.cases(core.d)]
site.d <- site.d[site %in% core.d$site]

#test fit with 3 plots.
#sub.site <- unique(core.d$site)[1:3]
#core.d <- core.d[site %in% sub.site]
#site.d <- site.d[site %in% sub.site]

#assign core level observations
   y <- core.d$relEM
c.x1 <- core.d$pH
core_site <- core.d$site           #this links core level observations to appropriate sites
core_site <- droplevels(core_site) #important to drop levels in subsetting data and losing sites. You will mess up the indexing otherwise.

#assign site level observations
s.x1 <- site.d$mat30
s.x2 <- site.d$map30
s.x3 <- site.d$doy
site_site <- site.d$site            #this indexes the site level variables and random effect of site.
site_site <- droplevels(site_site)

#transform [0,1] to (0,1) a la Cribari-Neto & Zeileis 2010
y <- (y * (length(y) - 1) + 0.5) / length(y)


jags.model = "
model{
# priors - a are core level, b are site level.
a0 ~ dnorm(0, .001)
a1 ~ dnorm(0, .001)
b1 ~ dnorm(0, .001)
b2 ~ dnorm(0, .001)
b3 ~ dnorm(0, .001)
tau <- pow(sigma,-2)
sigma ~ dunif(0, 100)
tau.site <- pow(sigma.site, -2)
sigma.site ~ dunif(0, 100)

#core level stuff - y is beta dsitributed.
for (i in 1:N){
      y[i] ~ dbeta(p[i], q[i])
      p[i] <- mu[i] * tau
      q[i] <- (1 - mu[i]) * tau
      #logit(mu[i]) <- a0 + site_effect[core_site[i]]
      logit(mu[i]) <- a1*c.x1[i] + site_effect[core_site[i]]
  }

#site level predictors
for(s in 1:S){
    #site_effect[s] <- random[site_site[s]] + b1*s.x1[site_site[s]]
    site_effect[s] <- b1*s.x1[site_site[s]] + b2*s.x2[site_site[s]] + b3*s.x3[site_site[s]] + random[site_site[s]]
    random[site_site[s]] ~ dnorm(0, tau.site) #distrubtion of random site effects.
  }

#forecast loop
#for (i in 1:N.new){
#y.new[i] ~ dbeta(p.new[i], q.new[i])
#p.new[i] <- mu.new[i] * tau
#q.new[i] <- (1 - mu.new[i]) * tau
#logit(mu.new[i]) <- a0 + a1*x1.new[i]
#}

} #end model
"

#list of data
jd <- list(y=y, N=length(core_site), core_site=core_site, c.x1=c.x1,                      #core level stuff
                S=length(site_site), site_site=site_site, s.x1=s.x1, s.x2=s.x2, s.x3=s.x3 #site level stuff
           )

#inital values chosen based on an LME model with same structure as JAGS model.
pre <- nlme::lme(boot::inv.logit(y) ~ pH + mat30 + map30 + doy, random = ~1|site, data = core.d)
inits <- list()
inits[[1]] <- list(a0 =  0.29, a1 = 0.01, b1 = -0.001, b2 = 0.000115, b3 = 0.00058)
inits[[2]] <- lapply(inits[[1]],'*',1.05)
inits[[3]] <- lapply(inits[[1]],'*',0.95)

#fit JAGS model. I can't get these chains to mix/converge when I add covariates.
jags.out <- autorun.jags(jags.model,
                     data=jd,
                     adapt = 1000,
                     #burnin = 20000,
                     #sample = 2000,
                     n.chains=3,
                     inits=inits,
                     method = 'rjparallel',
                     monitor=c('a0','a1', 'b1', 'b2','b3','random'))
#jags.out <- extend.jags(jags.out, sample = 5000, burnin = 10) #optionally extend the simulation if you haven't hit convergence.
summary(jags.out)
plot(jags.out)

saveRDS(jags.out,jags.out.path)
data.out <- list(site.d,core.d)
saveRDS(data.out,data.out.path)