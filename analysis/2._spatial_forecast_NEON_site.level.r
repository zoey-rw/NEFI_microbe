#Making a spatial forecast based on the prior data to NEON sites at the site level.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')


#load model results.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)

##### Quantify parameter covariance ####
pairs(mod$jags_model$mcmc[[1]], cex = 0.2)
ncol(mod$jags_model$mcmc[[1]])

z <- cov(mod$jags_model$mcmc[[1]])

#### Get NEON site level means and uncertainties ####

#### Setup data object ####
#1. needs predictor data and uncertainties.
#2. needs priors.