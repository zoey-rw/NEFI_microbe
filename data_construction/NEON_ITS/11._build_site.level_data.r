#Getting site-level NEON data for site-level prediction.
#clear environment, source paths.
rm(list = ls())
source('paths.r')
source('NEFI_functions/hierarch_core.means_JAGS.r')
source('NEFI_functions/hierarch_plot.means_JAGS.r')
source('NEFI_functions/crib_fun.r')

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
#to.ag <- to.ag[-207,]
relEM <- hierarch_plot.means_JAGS(x_mu = to.ag$relEM, plot_site = to.ag$siteID) #Can't figure out whats wrong with this.
test <- aggregate(basal_live ~ siteID, data = to.ag, FUN = mean)
test.sd <- aggregate(basal_live ~ siteID, data = to.ag, FUN = sd)
basal <- hierarch_plot.means_JAGS(x_mu = to.ag$basal_live, plot_site = to.ag$siteID)

#Okay- regular aggregate and bayesian aggregate giving different answers...
#But they give the same for percent C above...
z <- merge(basal$site.table, test, by = 'siteID')
plot(Mean ~ basal_live, data = z)
