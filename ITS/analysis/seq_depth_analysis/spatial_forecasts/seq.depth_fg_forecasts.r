#Making a spatial forecast for functional groups at multiple sequence depths.
#downstream this will log transform map, which prevents this from generalizing beyond the dirichlet example.
#This script depends on the following packages: DirichletReg.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')
source('NEFI_functions/ddirch_forecast.r')

#set output path.----
output.path <- ted_fg_seq.depth_ddirch_foreacsts.path

#load model results.----
mod_all <- readRDS(ted_fg_ddirch_fit_seqdepth.path)
#we need to loop the below code over each model in the "mod_all" object. Each model is some level of sequence depth.
super.list <- list()
#for(kk in 1:length(mod_all)){mod <- mod_all[[kk]]}

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

#Get forecasts from ddirch_forecast.----
core.fit <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID)
plot.fit <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID)
site.fit <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID)

#store output as a list and save.----
output <- list(core.fit,plot.fit,site.fit,core.preds,plot.preds,site.preds,core.sd,plot.sd,site.sd)
names(output) <- c('core.fit','plot.fit','site.fit',
                   'core.preds','plot.preds','site.preds',
                   'core.sd','plot.sd','site.sd')
saveRDS(output, output.path)
