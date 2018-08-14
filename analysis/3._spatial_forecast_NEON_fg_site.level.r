#Making a spatial forecast based on the prior data to NEON sites at the site level.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')

#load model results.
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)

#get site level covariates
site_covs <- readRDS(NEON_site_covs.path)
glob_covs <- readRDS(NEON_glob_covs.path)

n.samp <- 5000

for(i in 1:n.samp){
  mod <- mod[[1]]$jags_model
  preds <- as.character(mod[[1]]$species_parameter_output$other$predictor)
  covs <- site_covs[,grep(preds[3], colnames(site_covs))]
}


####1. sample rows of each model's parameters. This accounts for covariance.####

mod[[1]]$species_parameter_output

####2. Make predictions on the linear scale.####

####3. Convert predictions to the [0,1] scale.####

####4. repeat 2+3, setting parameter uncertainty to zero then covariate uncertainty to zero.####

####5. save all output for downstream plotting and variance decomposition. ####


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

