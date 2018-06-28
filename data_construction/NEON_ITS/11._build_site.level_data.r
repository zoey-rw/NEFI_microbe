#Getting site-level NEON data for site-level prediction.
#clear environment, source paths.
rm(list = ls())
source('paths.r')
source('NEFI_functions/hierarch_core.means_JAGS.r')
source('NEFI_functions/hierarch_plot.means_JAGS.r')

#load core level data.
core.dat <- readRDS(core.table.path)
core.dat <- core.dat[,c('siteID','plotID','sampleID','soilMoisture','soilInWaterpH')]

moisture <- hierarch_core.means_JAGS(x_mu = core.dat$soilMoisture , core_plot = core.dat$plotID)
      pH <- hierarch_core.means_JAGS(x_mu = core.dat$soilInWaterpH, core_plot = core.dat$plotID)

#load plot level data.
plot.dat <- readRDS(plot.table.path)
organicCpercent <- hierarch_plot.means_JAGS(x_mu = plot.dat$organicCPercent, plot_site = plot.dat$siteID)
        CNratio <- hierarch_plot.means_JAGS(x_mu = plot.dat$CNratio        , plot_site = plot.dat$siteID)
