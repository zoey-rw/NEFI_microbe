#Making a spatial forecast based on the prior data to NEON sites at the site level using dirichlet models.
#downstream this will log transform map, which prevents this from generalizing beyond the dirichlet example.
#This script depends on the following packages: DirichletReg.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')
source('NEFI_functions/ddirch_forecast_site.level.r')

#load model results.
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)

#get site and global level covariates
site_covs <- readRDS(NEON_site_covs.path)
site_sds  <- site_covs[,grep('_sd',colnames(site_covs))]
glob_covs <- readRDS(NEON_glob_covs.path)

#specify number of times to sample from parameters and covariates.
n.samp <- 1000

test <- ddirch_forecast_site.level(mod = mod,site_covs = site_covs,site_sds = site_sds,glob.covs = glob_covs)


#Here is where we gonna loop over models.
#this loop now runs in the function above.
prediction.output <- list()
for(i in 1:length(mod)){
  #grab the model out of list
  c.mod <- mod[[i]]
  j.mod <- c.mod$jags_model #separate out the jags model.

  #### organize your covariates and covariate standard deviations. ####
  preds <- as.character(c.mod$species_parameter_output$other$predictor)
  #grab covariates that were predictors.
  covs <- site_covs[,colnames(site_covs) %in% preds]
  #add intercept.
  covs <- cbind(rep(1,nrow(covs)), covs)
  colnames(covs)[1] <- 'intercept'
  #grab uncertainties in sd, if present.
  cov.sd <- precision_matrix_match(covs, site_sds)
  #re-order to match predictor order. Conditional otherwise it breaks the intercept only case.
  if(ncol(covs) > 1){
    covs   <-   data.frame(covs[,preds])
    cov.sd <- data.frame(cov.sd[,preds])
  }

  #fill in NAs from global level predictor means and sds.
  for(j in 1:ncol(covs)){
    cov.name <- colnames(covs)[j]
      covs[,j][is.na(  covs[,j])] <- glob_covs[glob_covs$predictor == cov.name,]$Mean
    cov.sd[,j][is.na(cov.sd[,j])] <- glob_covs[glob_covs$predictor == cov.name,]$SD
  }

  #### Sample parameter and covariate space, make predictions. ####
  pred.out <- list()
  cred.out <- list()
  for(j in 1:n.samp){
    #### Sample parameters from mcmc output. ####
    mcmc <- do.call(rbind,j.mod$mcmc)
    mcmc.sample <- mcmc[sample(nrow(mcmc),1),]
    #grab x.m values, convert to matrix.
    x.m <- mcmc.sample[grep("^x\\.m\\[", names(mcmc.sample))]
    x.m <- matrix(x.m, nrow = ncol(covs), ncol = length(x.m)/ncol(covs))

    #### Sample from covariate distributions. ####
    now.cov <- matrix(NA, ncol = ncol(covs), nrow = nrow(covs))
    for(k in 1:ncol(covs)){now.cov[,k] <- rnorm(nrow(covs),covs[,k], cov.sd[,k])}
    colnames(now.cov) <- preds
    now.cov <- data.frame(now.cov)
    #log transform map values if this is one of your covariates, put covariates back in matrix form.
    #anti-logit relEM, multiply by 100 if this is one of your covariates.
    if('relEM' %in% colnames(now.cov)){now.cov$relEM <- boot::inv.logit(now.cov$relEM) * 100}
    if('map'   %in% colnames(now.cov)){now.cov$map   <- log(now.cov$map)}
    now.cov <- as.matrix(now.cov)

    
    #### Combine covariates and parameters to make a prediction. ####
    pred.x.m <- matrix(NA, ncol=ncol(x.m), nrow = nrow(covs))
    for(k in 1:ncol(x.m)){pred.x.m[,k] <- exp(now.cov %*% x.m[,k])}
    #get mean prediction and then draw from dirichlet distribution.
    cred.out[[j]] <- pred.x.m / rowSums(pred.x.m)
    pred.out[[j]] <- DirichletReg::rdirichlet(nrow(pred.x.m),pred.x.m)
  }
  
  #### Summarize prediction mean and confidence intervals. ####
  pred.mean       <- apply(simplify2array(cred.out), 1:2, mean)
  pred.ci.0.025   <- apply(simplify2array(cred.out), 1:2, quantile, probs = c(0.025))
  pred.ci.0.975   <- apply(simplify2array(cred.out), 1:2, quantile, probs = c(0.975))
  pred.pi.0.025   <- apply(simplify2array(pred.out), 1:2, quantile, probs = c(0.025))
  pred.pi.0.975   <- apply(simplify2array(pred.out), 1:2, quantile, probs = c(0.975))
  
  #batch into list, column names are group names, row names are siteIDs.
  output <- list(pred.mean, pred.ci.0.025, pred.ci.0.975, pred.pi.0.025, pred.pi.0.975)
  for(k in 1:length(output)){
    colnames(output[[k]]) <- names(c.mod$species_parameter_output)
    rownames(output[[k]]) <- site_covs$siteID
  }
  names(output) <- c('mean','ci_0.025','ci_0.975','pi_0.025','pi_0.975')
  prediction.output[[i]] <- output
}#end fit loop.
names(prediction.output) <- names(mod)

####4. repeat 2+3, setting parameter uncertainty to zero then covariate uncertainty to zero.####

####5. save all output for downstream plotting and variance decomposition. ####
