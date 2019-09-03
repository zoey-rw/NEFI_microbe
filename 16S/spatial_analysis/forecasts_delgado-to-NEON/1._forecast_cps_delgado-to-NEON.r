
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/precision_matrix_match.r')
source('NEFI_functions/ddirch_forecast_noLogMap.r')

# source ddirch_forecast
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_forecast.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output path.----
# output.path <- NEON_cps_fcast_ddirch_16S.path
output.path <- paste0(scc_gen_16S_dir,"/NEON_forecasts/delgado-to-NEON_cps_forecast_ddirch.rds")

#load prior model results.----
all.mod <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_delgado/dir_delgado_8-30-19.rds"))

#get core-level covariate means and sd.----
dat <- readRDS(hierarch_filled_data.path)
dat <- lapply(dat, function(x) x[!(names(x) %in% c("pH", "conifer"))])
dat <- lapply(dat, function(x) setnames(x, old = "pH_water", new = "pH", skip_absent = TRUE))

core_mu <- dat$core.core.mu
plot_mu <- dat$plot.plot.mu
site_mu <- dat$site.site.mu
#site_mu$map <- log(site_mu$map) # the function already takes the log of map

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
fcast.output <- list()

cat('Making forecasts...\n')
for(i in 1:length(all.mod)){
 # for(i in 1:5){ # phylum only
 # for(i in 6:length(all.mod)){ # functional group only
    tic()
  mod <- all.mod[[i]]
  # core.fit <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID, n.samp = 1000)
  # plot.fit <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID  , n.samp = 1000)
  # site.fit <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID  , n.samp = 1000)
  core.fit <- ddirch_forecast_noLogMap(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID, n.samp = 1000)
  plot.fit <- ddirch_forecast_noLogMap(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID  , n.samp = 1000)
  site.fit <- ddirch_forecast_noLogMap(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID  , n.samp = 1000)
  
  #store output as a list and save.----
  output <- list(core.fit,plot.fit,site.fit,core.preds,plot.preds,site.preds,core.sd,plot.sd,site.sd)
  names(output) <- c('core.fit','plot.fit','site.fit',
                     'core.preds','plot.preds','site.preds',
                     'core.sd','plot.sd','site.sd')
  fcast.output[[i]] <- output
  cat(paste0(i,' of ',length(all.mod),' forecasts complete. '))
  toc()
}
cat('All forecasts complete.')
toc()


#Save output.----
names(fcast.output) <- names(all.mod)
saveRDS(fcast.output, output.path)
