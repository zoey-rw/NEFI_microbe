#Making a spatial forecast based on the prior data to NEON sites at the site level.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')


#load model results.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)

##### Quantify parameter means and covariance ####
#collapse 3 parameter chains to 1.
parms.all <- rbind(mod$jags_model$mcmc[[1]],mod$jags_model$mcmc[[2]],mod$jags_model$mcmc[[3]])
#calculate covariance
p.cov <- cov(parms.all)
#get parameter means and sd as vectors.
p.mu <- colMeans(parms.all)
p.sd <- apply(parms.all, 2, sd)

#Get a parameter draw from parameter distributions that accounts for covariance.
p.draw <- MASS::mvrnorm(p.mu,p.cov)


#Get NEON predictor matrix. Rows are sites, columns are x'values.
neon.t1 <- readRDS()
neon_core <- readRDS(core.table.path) #pH, moisture
neon_plot <- readRDS(plot.table.path) #%C, C:N
neon_site <- readRDS(site.table.path) #mat, map, mat_CV, map_CV, mdr, aridity, FOREST, relEM, CONIFER.

#aggregate cores to plots.
pH.plot_mu <- aggregate(soilInWaterpH ~ plotID, data = neon_core, FUN = mean)
pH.plot_mu$sd <- aggregate(soilInWaterpH ~ plotID, data = neon_core, FUN = sd)[,2]
pH.plot_mu$siteID <- substring(pH.plot_mu$plotID,1,4)
pH.site_mu <- aggregate(soilInWaterpH ~ siteID, data = pH.plot_mu, FUN = mean)
pH.site_mu$pH_sd <- aggregate(soilInWaterpH ~ siteID, data = pH.plot_mu, FUN = sd)[,2]
pH.plot_mu[is.na(pH.plot_mu$sd),]$sd <- mean(pH.plot_mu$sd, na.rm = T)
siteID <- pH.plot_mu$siteID

#draw from plot-level uncertainty
n.sim <- 1000
out.list <- list()
for(i in 1:n.sim){
  plot.mu <- vector(length = nrow(pH.plot_mu))
  for(k in 1:nrow(pH.plot_mu)){
    plot.mu[k] <- rnorm(1,pH.plot_mu$soilInWaterpH[k],pH.plot_mu$sd[k])  
  }
  data <- data.frame(cbind(plot.mu, siteID))
  data$plot.mu <- as.numeric(data$plot.mu)
  mu   <- aggregate(data$plot.mu ~ data$siteID, FUN = mean)
  out.list[[i]] <- mu[,2]
}
test <- do.call('rbind',out.list)
colnames(test) <- mu[,1]

#Just setup a JAGS model that generates site level values of everything.
jags.model="
model {
  #take means hierarchically.
  for(i in 1:n.core){
    plot_mu[] ~ dnorm()
  }
  #fill in any missing X values at core level.
  for(m in 1:N.core.preds){
    for(i in 1:N.site){ site.core[i,m] ~ dnorm(core.global[m], core.site.tau[m])}           ## estimate global mean
    for(j in 1:N.plot){ plot.core[j,m] ~ dnorm(site.core[plot_site[j],m],core.plot.tau[m])} ## estimate site means
    for(k in 1:N.core){   core_mu[k,m] ~ dnorm(plot.core[core_plot[k],m],core.core.tau[m])} ## estimate plot means
  }
}
"

