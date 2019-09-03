#Making a spatial forecast based on the prior data to NEON sites at core, plot and site levels.
#forecasting functional groups.
#downstream this will log transform map, which prevents this from generalizing beyond the dirichlet example.
#This script depends on the following packages: DirichletReg.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/precision_matrix_match.r')
source('NEFI_functions/tic_toc.r')
#source('NEFI_functions/ddirch_forecast.r')
library(data.table)

library(RCurl)
# source paths.r from colins github
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/paths.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

compl_cases <- F
if (compl_cases== T){
  output.path <- paste0(pecan_gen_16S_dir,"NEON_forecast_data/NEON_fcast_fg_comp_cases.rds")
} else output.path <- NEON_cps_fcast_fg_16S.path

#load prior model fits----
fg <- readRDS(prior_16S_all.fg.groups_JAGSfits.path)
#phylo <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/JAGS_output/bahram_16S.prior_phylo_new_test.rds")
## combine two lists of models 
#all_mods <- do.call(c, list(phylo,fg))

all_fcasts <- list()
#all_fcasts <- readRDS(NEON_cps_fcast_fg_16S.path)
for (p in 1:length(fg)) {
     #12:length(fg)) {
mod <- fg[[p]]
mod <- mod$no.nutr.preds
#mod <- mod$all.preds #just the selected covariates

#get core-level covariate means and sd.----
if (compl_cases==T) { dat <- readRDS(missing_data_removed_16S.path)
}else dat <- readRDS(hierarch_filled_16S.path) # using ITS data right now.
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

# remove micronutrients.
lst <- list(core.preds, core.sd, plot.preds, plot.sd, site.preds, site.sd)
names(lst) <- c("core.preds", "core.sd", "plot.preds", "plot.sd", "site.preds", "site.sd")
list2env(
  lapply(lst, function(x) x[!(names(x) %in% c("P", "K", "Mg", "Ca"))]), 
         envir=.GlobalEnv)


#Get forecasts from ddirch_forecast.----
tic()
core.fit <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID)
plot.fit <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID)
site.fit <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID)
toc()
cat(paste("Forecast created for model",p,"\n"))

#store output as a list and save.----
output <- list(core.fit,plot.fit,site.fit,core.preds,plot.preds,site.preds,core.sd,plot.sd,site.sd)
names(output) <- c('core.fit','plot.fit','site.fit',
                   'core.preds','plot.preds','site.preds',
                   'core.sd','plot.sd','site.sd')
all_fcasts[[p]] <- output
}

saveRDS(all_fcasts, output.path)
