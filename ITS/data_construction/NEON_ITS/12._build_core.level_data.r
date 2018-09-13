#build core-level observation data set.
#clear environment, source paths.
rm(list=ls())
library(runjags)
source('paths.r')
source('NEFI_functions/crib_fun.r')
source('NEFI_functions/pC_uncertainty_neon.r')
source('NEFI_functions/cn_uncertainty_neon.r')

#1. observations made at the core level.----
#samples with ITS data, along with DNA identifiers. Only take samples measured during peak greeness.
dp1.10801 <- readRDS(dp1.10108.00_output.path)
dp1.10801$geneticSampleID <- as.character(dp1.10801$geneticSampleID)
dp1.10801 <- dp1.10801[!(dp1.10801$geneticSampleID == ''),] #many samples missing geneticSampleID. throw them out.
dp1.10801 <- dp1.10801[order(dp1.10801$geneticSampleID),]

#moisture and pH data.
dp1.10086 <- readRDS(dp1.10086.00_output.path)
dp1.10086$site_date <- paste0(dp1.10086$siteID,'-',dp1.10086$dateID)
dp1.10086$geneticSampleID <- as.character(dp1.10086$geneticSampleID)

#soil C and N data.
dp1.10078 <- readRDS(dp1.10078.00_output.path)
dp1.10078$site_date_plot <- paste0(dp1.10078$siteID,'-',dp1.10078$dateID,'-',dp1.10078$plotID)

#merge these 3 sets of observations together.
to_merge <- dp1.10086[,!c(colnames(dp1.10086) %in% colnames(dp1.10801))]
to_merge$geneticSampleID <- dp1.10086$geneticSampleID
merged <- merge(dp1.10801,to_merge)
to_merge <- dp1.10078[,!c(colnames(dp1.10078) %in% colnames(merged))]
to_merge$sampleID <- dp1.10078$sampleID
to_merge <- to_merge[!(duplicated(to_merge$sampleID)),] #get rid of analytical replicates.
merged <- merge(merged,to_merge, by = 'sampleID', all.x=T)
merged$year <- substring(merged$dateID,1,4)
merged <- merged[!is.na(merged$siteID),]

#Peak Greenness in 2014. 551 observations. 13 sites.
merged <- merged[merged$sampleTiming == 'peakGreenness' & merged$year == '2014',]

#finalize columns for core.level.
core.level <- merged[,c('sampleID','geneticSampleID','dnaSampleID','siteID','plotID','dateID','collectDate','horizon','elevation','soilMoisture','soilInWaterpH','organicCPercent','CNratio')]
colnames(core.level)[(ncol(core.level) - 2) : ncol(core.level)] <- c('pH','pC','cn')
core.level$pC_sd <- pC_uncertainty_neon(core.level$pC)
core.level$cn_sd <- cn_uncertainty_neon(core.level$cn)

#2. plot level observations.----
#relative abundance EM trees
dp.10098.plot <- readRDS(dp1.10098.00_plot.level.path)


to.ag <- dp.10098.plot[,c('plotID','relEM','basal_live')]
#gotta crib and logit transform to scale up.
to.ag$relEM <- boot::logit(crib_fun(to.ag$relEM))
plot.level <- aggregate(.~plotID, to.ag, mean, na.rm=T)
plot_sd <- aggregate(.~plotID, to.ag, sd, na.rm = T)
plot.level$siteID <- substring(plot.level$plotID,1,4)

#3. site level observations.----
site.level <- readRDS(site_level_data.path)



