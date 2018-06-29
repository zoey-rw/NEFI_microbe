#getting hierarchical means c/p/s
rm(list=ls())
library(runjags)

#### pseudo-data generation ####
#4 sites with their own means and standard deviations.
site_mu <- c(20,30,19,12)
site_sd <- c(2,4,2,1)
siteID <- c(LETTERS[1:4])

#specify number of plots per site and number of cores per plot.
n.plot <- 20
n.core <- 30
  
plots <- list()
for(i in 1:length(site_mu)){
  plot_mu <- rnorm(n.plot,site_mu[i],site_sd[i])
  plots[[i]]  <- plot_mu
}
plot_mu <- unlist(plots)
plotID <- list()
for(i in 1:length(siteID)){
    plots <- vector(length = n.plot)
    for(k in 1:n.plot){
      plots[k] <- paste0(siteID[i],'_',k)
    }
  plotID[[i]] <- plots
}
plotID <- unlist(plotID)
siteID <- substring(plotID,1,1)
dat <- data.frame(siteID,plotID,plot_mu)

core_sd <- 1
cores <- list()
for(i in 1:length(plot_mu)){
  core_mu <- rnorm(n.core,plot_mu[i],core_sd)
  cores[[i]] <- core_mu
}
core_mu <- unlist(cores)
coreID <- list()
for(i in 1:length(plotID)){
  cores <- vector(length=n.core)
  for(k in 1:n.core){
    cores[k] <- paste0(plotID[i],"_",k)
  }
  coreID[[i]] <- cores
}
coreID <- unlist(coreID)
plotID <- sub('_([^_]*)$', '', coreID)
siteID <- substring(plotID,1,1)
dat <- data.frame(core_mu,coreID,plotID,siteID)

#### specify JAGS model ####
jags.model = "
model {
  #get plot level means
  for(i in 1:N.core){
    core_mu[i] ~ dnorm(plot_mu[core_plot[i]], tau.plot[core_plot[i]])
  }

  #get site level means.
  for(i in 1:N.plot) {
    plot_mu[i] ~ dnorm(site_mu[plot_site[i]], tau.site[plot_site[i]])
  }

  #get global mean and uncertainty.
  for(i in 1:N.site){
    site_mu[i] ~ dnorm(global_mu, tau.glob)
  }

  #Plot level priors
  for(i in 1:N.plot){
    sigma.plot[i] ~ dunif(0, 100)
      tau.plot[i] <- pow(sigma.plot[i], -2)
  }
  #Site level priors
  for(i in 1:N.site) {
    sigma.site[i] ~ dunif(0, 100)
      tau.site[i] <- pow(sigma.site[i], -2)
  }
  #Global level priors
   global_mu ~ dnorm(0,1E-3)
  sigma.glob ~ dunif(0,100)
    tau.glob <- pow(sigma.glob, -2)
    
}"

#### setup JAGS data object. ####
jags.data <- list(N.core = nrow(dat), N.plot = length(unique(dat$plotID)), N.site = length(unique(dat$siteID)),
                  core_mu = dat$core_mu, core_plot = droplevels(as.factor(dat$plotID)), plot_site = droplevels(as.factor(substring(unique(dat$plotID), 1,1)))
                  )

#### fit JAGS model ####
test <- run.jags(model = jags.model,
                  data = jags.data,
              n.chains = 3,
               monitor = c('site_mu','global_mu'),
                 adapt = 200,
                burnin = 1000,
                sample = 1000)
out <- summary(test)
