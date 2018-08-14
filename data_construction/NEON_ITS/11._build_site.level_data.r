#Getting site-level NEON data for site-level prediction.
#clear environment, source paths.
rm(list = ls())
library(runjags)
source('paths.r')
source('NEFI_functions/hierarch_core.means_JAGS.r')
source('NEFI_functions/hierarch_plot.means_JAGS.r')
source('NEFI_functions/crib_fun.r')

#specify output paths
site_output_path <- NEON_site_covs.path
glob_output.path <- NEON_glob_covs.path

#load core level data.
core.dat <- readRDS(core.table.path)
core.dat <- core.dat[,c('siteID','plotID','sampleID','soilMoisture','soilInWaterpH')]

moisture <- hierarch_core.means_JAGS(x_mu = core.dat$soilMoisture , core_plot = core.dat$plotID)
      pH <- hierarch_core.means_JAGS(x_mu = core.dat$soilInWaterpH, core_plot = core.dat$plotID)

#load plot level data.
plot.dat <- readRDS(plot.table.path)
organicCpercent <- hierarch_plot.means_JAGS(x_mu = plot.dat$organicCPercent, plot_site = plot.dat$siteID)
        CNratio <- hierarch_plot.means_JAGS(x_mu = plot.dat$CNratio        , plot_site = plot.dat$siteID)


#Get relative abundance EM
dp.10098.plot <- readRDS(dp1.10098.00_plot.level.path)
to.ag <- dp.10098.plot[,c('plotID','siteID','relEM','basal_live')]
#gotta crib and logit transform to scale up.
to.ag$relEM <- boot::logit(crib_fun(to.ag$relEM))
relEM <- hierarch_plot.means_JAGS(x_mu = to.ag$relEM, plot_site = to.ag$siteID, gamma_variance = T)
basal <- hierarch_plot.means_JAGS(x_mu = to.ag$basal_live, plot_site = to.ag$siteID)

#Get wheter its a forest or not.Colin just got this from the site descriptions, data is below.
#Consider dropping STER, its agriculture.
siteID <- c('DSNY','HARV','BART','JERC','ORNL','SCBI','TALL','UNDE','CPER','STER','WOOD')
 forest <- c(0,1,1,0,1,1,1,1,0,0,0)
conifer <- c(1,1,1,0,1,0,1,1,0,0,0)
forest.data <- data.frame(siteID,forest,conifer)

#get site level data products
site.dat <- readRDS(site_level_data.path)

#begin merging site level.
out <- merge(site.dat, forest.data, all.x=T)
out <- merge(out, relEM$site.table[,c('Mean','SD','siteID')], all.x=T)
names(out)[names(out) == 'Mean'] <- 'relEM'
names(out)[names(out) == 'SD'  ] <- 'relEM_sd'
out <- merge(out, basal$site.table[,c('Mean','SD','siteID')], all.x=T)
names(out)[names(out) == 'Mean'] <- 'basal'
names(out)[names(out) == 'SD'  ] <- 'basal_sd'
out <- merge(out, organicCpercent$site.table[,c('Mean','SD','siteID')], all.x=T)
names(out)[names(out) == 'Mean'] <- 'pC'
names(out)[names(out) == 'SD'  ] <- 'pC_sd'
out <- merge(out, CNratio$site.table[,c('Mean','SD','siteID')], all.x=T)
names(out)[names(out) == 'Mean'] <- 'cn'
names(out)[names(out) == 'SD'  ] <- 'cn_sd'
out <- merge(out, moisture$site.table[,c('Mean','SD','siteID')], all.x=T)
names(out)[names(out) == 'Mean'] <- 'moisture'
names(out)[names(out) == 'SD'  ] <- 'moisture_sd'
out <- merge(out, pH$site.table[,c('Mean','SD','siteID')], all.x=T)
names(out)[names(out)=='Mean'] <- 'pH'
names(out)[names(out) == 'SD'  ] <- 'pH_sd'

#deal with grassland NAs for relEM and basal area.
out[out$forest == 0 & is.na(out$basal   ),]$basal    <- 0
out[out$forest == 0 & is.na(out$basal_sd),]$basal_sd <- 0.01
out[out$forest == 0 & is.na(out$relEM   ),]$relEM    <- -10
out[out$forest == 0 & is.na(out$relEM_sd),]$relEM_sd <- 0.01

if(out$forest == 0 & is.na(out$basal)){
  out$basal <- 0
}


#get global level uncertainties where relevant.
glob.out <-
rbind(             pH$glob.table[,c('Mean','SD')],
             moisture$glob.table[,c('Mean','SD')],
              CNratio$glob.table[,c('Mean','SD')],
      organicCpercent$glob.table[,c('Mean','SD')],
                basal$glob.table[,c('Mean','SD')],
                relEM$glob.table[,c('Mean','SD')])
glob.out$predictor <- c('pH','moisture','cn','pC','basal','relEM')


#save output
saveRDS(     out, site_output_path)
saveRDS(glob.out, glob_output.path)