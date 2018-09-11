#hierarchically aggregate multivariate composition data c/p/s using dirchlet.
#start 3 groups, 4 sites, plots w/in site, cores w/in plots.
#works ignoring plots within sites.
#Next add plots in between core and site level.
rm(list=ls())
library(runjags)

#1. generate pseudo data.----
#3 species relative abundances across 4 sites.
spp1 <- c(0.2,0.8,0.5,0.1)
spp2 <- c(0.3,0.1,0.4,0.1)
spp3 <- 1 - (spp1 + spp2)
y <- data.frame(spp1,spp2,spp3)
y <- as.matrix(y)
site_sd <- 0.02

#get some numbevr of pltos and cores within plots.
n.site <- nrow(y)
n.plot <- 15
n.core <- 10

#simulate from dirichlet distribution.
site.list <- list()
for(i in 1:n.site){
  plot.list <- list()
  for(j in 1:n.plot){
      cores <- DirichletReg::rdirichlet(n.core,y[i,])
      cores <- data.frame(cbind(rep(LETTERS[i],n.core),rep(LETTERS[j],n.core),LETTERS[1:n.core],cores))
      colnames(cores) <- c('siteID','plotID','coreID','spp1','spp2','spp3')
      plot.list[[j]] <- cores
  }
  site.list[[i]] <- do.call(rbind,plot.list)
}
dat <- do.call(rbind,site.list)
dat[,4:6] <- sapply(dat[4:6],as.character)
dat[,4:6] <- sapply(dat[4:6],as.numeric)
y <- dat[,4:6]

#core_plot, plot_site factors.
dat$plotID <- paste0(dat$siteID,'_',dat$plotID)
core_plot = dat$plotID
plot_site = dat[seq(1, nrow(dat), 3), ]$siteID


#2. specify JAGS model.----
jags.model = "
model {
#get plot level means.
for(i in 1:N.core){
  for(j in 1:N.spp){
    log(core.hat[i,j]) <- alpha + plot_mu[core_plot[i],j]
  }
  y[i,1:N.spp] ~ ddirch(core.hat[i,1:N.spp]) 
}

#get site level means.
for(i in 1:N.plot){
  for(j in 1:N.spp){
    #log(plot.hat[i,j]) <- alpha + site_mu[plot_site[i],j]
    log(plot.hat[i,j]) <- site_mu[plot_site[i],j]
  }
  plot_mu[i,1:N.spp] ~ ddirch(plot.hat[i,1:N.spp])
}

#prior
#for(i in 1:N.site){for(j in 1:N.spp){site_mu[i,j] ~ dnorm(0,1E-3)}}
alpha ~ dnorm(0, 1.0E-3) 
for(i in 1:N.site){
  site_mu[i,1] <- 0
  for(j in 2:N.spp){
    site_mu[i,j] ~ dnorm(0,1E-3)
    }
}

#map to something useful.
#for(i in 1:N.site){
#  final_mu[i,1] <- alpha
#    for(j in 2:N.spp){
#      final_mu[i,j] <- site_mu[i,j] + alpha
#    }
#  }
}"

#3. Setup JAGS data object.----
jd <- list(y=as.matrix(y), N.site = n.site, N.core = nrow(y), N.plot = length(unique(dat$plotID)), N.spp = ncol(y), 
           core_plot = as.factor(core_plot),
           plot_site = as.factor(plot_site))

#4. Run JAGS model.----
test <- run.jags(model = jags.model,
                 data = jd,
                 n.chains = 3,
                 monitor = c('site_mu'),
                 adapt = 200,
                 burnin = 2000,
                 sample = 1000)
out <- summary(test)
spp.sum <- matrix(out[,4], nrow = n.site, ncol = ncol(y))
spp.sum <- boot::inv.logit(spp.sum) / rowSums(boot::inv.logit(spp.sum))
spp.sum
