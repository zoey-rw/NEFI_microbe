#Fit a hierarchical core/plot/site jags model.
#clearn envirionment, source paths.
rm(list=ls())
source('paths.r')
#source('/home/caverill/NEFI_microbe/data_formatting/formatting_NEON_microbial/core_site_plot_aggregation_May.2018/0. aggregation paths.r')
#source some other functions
source('NEFI_functions/hierarchical_linear_dirlichet_jags.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/mean_center.r')

#load data.
obs.table <- readRDS( obs.table.path)
core.table <- readRDS(core.table.path)
plot.table <- readRDS(plot.table.path)
site.table <- readRDS(site.table.path)

##TESTING
#lets only look at complete cases for 2 sites right now.
obs.table <-  obs.table[obs.table$siteID %in% c('BART','HARV'),]
core.table <- core.table[core.table$geneticSampleID %in% obs.table$geneticSampleID,]
plot.table <- plot.table[plot.table$plotID %in% obs.table$plotID,]
site.table <- site.table[site.table$siteID %in% obs.table$siteID,]


#put everything in correct order. 
#Convert all factor vectors to character vectors - this line sucks and is a good argument for data.table/tidyverse
 obs.table[sapply( obs.table, is.factor)] <-  lapply( obs.table[sapply( obs.table, is.factor)], as.character)
core.table[sapply(core.table, is.factor)] <-  lapply(core.table[sapply(core.table, is.factor)], as.character)
plot.table[sapply(plot.table, is.factor)] <-  lapply(plot.table[sapply(plot.table, is.factor)], as.character)
site.table[sapply(site.table, is.factor)] <-  lapply(site.table[sapply(site.table, is.factor)], as.character)

#hierarchically order the rows as site/plot/core.
obs.table <-  obs.table[order( obs.table$siteID, obs.table$plotID, obs.table$geneticSampleID),]
core.table <- core.table[order(core.table$siteID,core.table$plotID,core.table$geneticSampleID),]
plot.table <- plot.table[order(plot.table$siteID,plot.table$plotID),]
site.table <- site.table[order(site.table$siteID),]
core_mu <- core.table
plot_mu <- plot.table
site_mu <- site.table

#setup indexing
core_plot <- as.factor(core.table$plotID)
core_site <- as.factor(core.table$siteID)
plot_plot <- as.factor(plot.table$plotID)
site_site <- as.factor(site.table$siteID)

#setup y, core, plot and site data frames, and associated sd tables where applicable.
y <- obs.table[,3:6]
#dropping soil moisture and soil pH for now because all NA values for two site-dates of interest.
#core <- core.table[,c('siteID','plotID','soilMoisture','soilInWaterpH','soilTemp')]
#colnames(core)[3:5] <- c('moisture','pH','temperature')
core <- core.table[,c('siteID','plotID','soilTemp')]
colnames(core)[3] <- c('temperature')
core_sd <- core.table[,c('siteID','plotID')]
plot <- plot.table[,c('siteID','plotID','organicCPercent','CNratio')]
plot_sd <- plot.table[,c('siteID','plotID','organicCPercent_sd','CNratio_sd')]
site <- site.table[,c('siteID','elevation','map','mat','n.dep')]
site_sd <- site.table[,c('siteID','map_sd','mat_sd')]

#drop plot-site indices from data frames. Those saved in index vectors.
core$intercept <- rep(1,nrow(core))
core <- core[,-c(1:2)]; core <- core[,c('intercept','temperature')]
plot <- plot[,-c(1:2)]
site <- site[,-c(1:2)]
for(i in 1:ncol(core)){
  core[,i] <- mean_center(core[,i])
}
for(i in 1:ncol(plot)){
  plot[,i] <- mean_center(plot[,i])
}
for(i in 1:ncol(site)){
  site[,i] <- mean_center(site[,i])
}

core_sd <- NA
plot_sd <- NA
site_sd <- NA


#transform y values.
y <- data.frame(lapply(y, crib_fun))
#reorder y columns so 'other' is first.
y <- data.frame(y$other, y[,1:3])

core_mu = core; plot_mu = plot; site_mu = site

test <- hierarchical_dirlichet_jags.r(y = y, 
                                      core_mu = core, plot_mu = plot, site_mu = site,
                                      core_sd = core_sd, plot_sd = plot_sd, site_sd = site_sd,
                                      core_plot = core_plot, core_site = core_site, plot_plot = plot_plot, site_site = site_site)
