#spatial forecast variance decomposition: fungal functional groups.
#clear envrionment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/ddirch_forecast.r')

#load model results.----
#mod 1 is data from maps.
#mod 2 is site-specific data, no maps.
#mod 3 is all covariates.
mod <- readRDS(ted_ITS.prior_fg_JAGSfit)
mod <- mod[[3]] #just the all predictor case.

#get core-level covariate means and sd.----
dat <- readRDS(hierarch_filled.path)
core_mu <- dat$core.core.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#merge together.
plot_mu$siteID <- NULL
core.preds <- merge(core_mu   , plot_mu)
core.preds <- merge(core.preds, site_mu)
core.preds$relEM <- NULL
names(core.preds)[names(core.preds)=="b.relEM"] <- "relEM"

#get core-level SD.
core_sd <- dat$core.core.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together.
plot_sd$siteID <- NULL
core.sd <- merge(core_sd   , plot_sd)
core.sd <- merge(core.sd, site_sd)
core.sd$relEM <- NULL
names(core.sd)[names(core.sd)=="b.relEM"] <- "relEM"


#get plot-level covariate means and sd.----
core_mu <- dat$core.plot.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#merge together, .
plot_mu$siteID <- NULL
plot.preds <- merge(core_mu,plot_mu)
plot.preds <- merge(plot.preds,site_mu)
plot.preds$relEM <- NULL
names(plot.preds)[names(plot.preds)=="b.relEM"] <- "relEM"

#get plot-level SD.
core_sd <- dat$core.plot.sd
plot_sd <- dat$plot.plot.sd
site_sd <- dat$site.site.sd
#merge together
plot_sd$siteID <- NULL
plot.sd <- merge(core_sd,plot_sd)
plot.sd <- merge(plot.sd,site_sd)
plot.sd$relEM <- NULL
names(plot.sd)[names(plot.sd)=='b.relEM'] <- "relEM"

#get site-level covariate means and sd.----
core_mu <- dat$core.site.mu
plot_mu <- dat$plot.site.mu
site_mu <- dat$site.site.mu
#merge together
site.preds <- merge(core_mu, plot_mu)
site.preds <- merge(site.preds,site_mu)
names(site.preds)[names(site.preds)=='b.relEM'] <- "relEM"

#get site-level SD.
core_sd <- dat$core.site.sd
plot_sd <- dat$plot.site.sd
site_sd <- dat$site.site.sd
#merge together.
site.sd <- merge(core_sd,plot_sd)
site.sd <- merge(site.sd,site_sd)
names(site.sd)[names(site.sd)=='b.relEM'] <- "relEM"

#test the forecast function. works.
ref <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID)
zero_covariate <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID, zero_covariate_uncertainty = T)
zero_parameter <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID, zero_parameter_uncertainty = T)
  zero_process <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID, zero_process_uncertainty = T)

#try with plot-level.
ref <- ddirch_forecast(mod=mod, cov_mu = plot.preds, cov_sd = plot.sd, names = plot.preds$plotID)
zero_covariate <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID, zero_covariate_uncertainty = T)
zero_parameter <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID, zero_parameter_uncertainty = T)
  zero_process <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID,   zero_process_uncertainty = T)

#try with site-level.
ref <- ddirch_forecast(mod=mod, cov_mu = site.preds, cov_sd = site.sd, names = site.preds$siteID)
zero_covariate <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID, zero_covariate_uncertainty = T)
zero_parameter <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID, zero_parameter_uncertainty = T)
  zero_process <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID,   zero_process_uncertainty = T)
  
  

#generally goes down, the biggest effect is process.
#Prob need more simulations to really lock down small differences due to parameter and process.
head(ref$pi_0.975)
head(zero_covariate$pi_0.975)
head(zero_parameter$pi_0.975)
head(zero_process$pi_0.975)

#testing function
#mod=mod;cov_mu=core.preds;cov_sd=core.sd;names=core.preds$sampleID;n.samp = 1000
#zero_parameter_uncertainty=F
#zero_covariate_uncertainty=T
