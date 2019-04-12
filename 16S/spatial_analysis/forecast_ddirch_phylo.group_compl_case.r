#Making a 16S spatial forecast based on the prior data to NEON sites at core, plot and site levels.
# All phylogenetic levels. All taxa that are present in more than 50% of prior sites.
#downstream this will log transform map, which prevents this from generalizing beyond the dirichlet example.
#This script depends on the following packages: DirichletReg.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
source('NEFI_functions/tic_toc.r')
source('NEFI_functions/precision_matrix_match.r')
#source('NEFI_functions/ddirch_forecast.r')

# source ddirch_forecast
library(RCurl)
script <- getURL("https://raw.githubusercontent.com/colinaverill/NEFI_microbe/master/NEFI_functions/ddirch_forecast.r", ssl.verifypeer = FALSE)
eval(parse(text = script))

#set output path.----
output.path <- paste0(pecan_gen_16S_dir,"/NEON_forecast_data/NEON_fcast_comp_cases.rds")

#load prior model results.----
#all.mod <- readRDS(paste0(scc_gen_16S_dir,'JAGS_output/bahram_16S.prior_phylo_new_test.rds'))
#all.mod <- readRDS(bahram_16S_prior_phylo.group_JAGSfits)
#all.mod <- readRDS(paste0(scc_gen_16S_dir,"/JAGS_output/prior_phylo_JAGSfit_phylum.rds"))
all.mod <- readRDS(paste0(scc_gen_16S_dir, "JAGS_output/prior_phylo_JAGSfit_fewer_taxa.rds"))
phylum.mod <- readRDS(paste0(scc_gen_16S_dir,"JAGS_output/prior_phylo_JAGSfit_phylum_fewer_taxa_more_burnin.rds"))
all.mod$phylum <- phylum.mod$phylum

#get core-level covariate means and sd.----
dat <- readRDS(missing_data_removed_16S.path)
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
phylo.output <- list()

cat('Making forecasts...\n')
for(i in 1:length(all.mod)){
  tic()
  mod <- all.mod[[i]]
  mod <- mod$no.nutr.preds
  core.fit <- ddirch_forecast(mod=mod, cov_mu=core.preds, cov_sd=core.sd, names=core.preds$sampleID, n.samp = 1000)
  plot.fit <- ddirch_forecast(mod=mod, cov_mu=plot.preds, cov_sd=plot.sd, names=plot.preds$plotID  , n.samp = 1000)
  site.fit <- ddirch_forecast(mod=mod, cov_mu=site.preds, cov_sd=site.sd, names=site.preds$siteID  , n.samp = 1000)
  
  #store output as a list and save.----
  output <- list(core.fit,plot.fit,site.fit,core.preds,plot.preds,site.preds,core.sd,plot.sd,site.sd)
  names(output) <- c('core.fit','plot.fit','site.fit',
                     'core.preds','plot.preds','site.preds',
                     'core.sd','plot.sd','site.sd')
  phylo.output[[i]] <- output
  cat(paste0(i,' of ',length(all.mod),' forecasts complete. '))
  toc()
}
cat('All forecasts complete.')
toc()

#Save output.----
names(phylo.output) <- names(all.mod)
saveRDS(phylo.output, output.path)

# phylum: saveRDS(phylo.output, paste0(pecan_gen_16S_dir,"/NEON_forecast_data/NEON_fcast_phylum.rds"))

