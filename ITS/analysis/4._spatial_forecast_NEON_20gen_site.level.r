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
mod <- readRDS(ted_ITS.prior_20gen_JAGSfit)

#get site and global level covariates
site_covs <- readRDS(NEON_site_covs.path)
site_sds  <- site_covs[,grep('_sd',colnames(site_covs))]
glob_covs <- readRDS(NEON_glob_covs.path)

#specify number of times to sample from parameters and covariates.
n.samp <- 3000

fcast <- ddirch_forecast_site.level(mod = mod,site_covs = site_covs,site_sds = site_sds,glob.covs = glob_covs)

#save output.
saveRDS(fcast,NEON_site_fcast_20gen.path)
