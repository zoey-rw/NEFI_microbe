#Core level aggregation.
#Getting core, plot, site and global means and sd's for all observations made at the core level.
#This is currently soil pH, soil %C and soil C:N.
#clear environment, source paths.
rm(list=ls())
library(runjags)
library(dplyr)
source('paths.r')
source('NEFI_functions/pC_uncertainty_neon.r')
source('NEFI_functions/cn_uncertainty_neon.r')
source('NEFI_functions/hierarch_core.means_JAGS.r')
source('NEFI_functions/convert_K.r')

# load pH values for conversion
pH_conversion <- readRDS(pH_conversion.path)

#load data, subset and order.----
#load samples with ITS *or* 16S data, along with DNA identifiers.
dp1.10801 <- readRDS(dp1.10108.001_output.path)
dp1.10801 <- dp1.10801$core_16S_ITS
dp1.10801$geneticSampleID <- as.character(dp1.10801$geneticSampleID)
dp1.10801 <- dp1.10801[!(dp1.10801$geneticSampleID == ''),] #many samples missing geneticSampleID. throw them out.
dp1.10801 <- dp1.10801[order(dp1.10801$geneticSampleID),]

#moisture and pH data.
dp1.10086 <- readRDS(dp1.10086.001_output.path)
dp1.10086$site_date <- paste0(dp1.10086$siteID,'-',dp1.10086$dateID)
dp1.10086$geneticSampleID <- as.character(dp1.10086$geneticSampleID)

#soil C and N data.
dp1.10078 <- readRDS(dp1.10078.001_output.path)
# remove analytical replicates
dp1.10078 <- dp1.10078[-which(dp1.10078$remarks %in% c("Duplicate samples", "Replicate samples", "Replicate samples.")),]
dp1.10078 <- dp1.10078[-which(dp1.10078$analyticalRepNumber == 2),]

# fix chem data with N and C split up
dup <- dp1.10078[duplicated(dp1.10078$sampleID),]
dup <- dp1.10078[dp1.10078$sampleID %in% dup$sampleID,] # we want all duplicated rows
N.only <- dup[dup$acidTreatment == "N",]
C.only <- dup[dup$acidTreatment == "Y",]
C.only$nitrogenPercent <- NULL
N.only.merge <- N.only[,colnames(N.only) %in% c("sampleID","nitrogenPercent")]
collapsed.C.N <- merge(C.only, N.only.merge)

# merge back in with rest of C/N data
dp1.10078 <- dp1.10078[!(dp1.10078$sampleID %in% collapsed.C.N$sampleID),]
dp1.10078$CNvals_merged <- FALSE
collapsed.C.N$CNvals_merged <- TRUE # so we know for later that these values are different
collapsed.C.N <- collapsed.C.N[,colnames(dp1.10078)] # reorder column names
collapsed.C.N$CNratio <- round(collapsed.C.N$organicCPercent/collapsed.C.N$nitrogenPercent,1)
dp1.10078 <- rbind(dp1.10078, collapsed.C.N) # merge back together


#merge phys and chem, then DNA observations together.----
dp1.10078.merge <- dp1.10078[,!c(colnames(dp1.10078) %in% colnames(dp1.10086))]
keep_cols <- c("sampleID", "namedLocation", "siteID", "plotID", "collectDate", "dateID")
dp1.10078.merge[,keep_cols] <- dp1.10078[,keep_cols]
merged <- merge(dp1.10086,dp1.10078.merge, all.x=T) # physical data + chemistry data

# merge phys/chem with DNA observations
dp1.10801.merge <- dp1.10801[,!c(colnames(dp1.10801) %in% colnames(merged))]
keep_cols <- c("geneticSampleID", "namedLocation", "siteID", "plotID", "collectDate", "dateID")
dp1.10801.merge[,keep_cols] <- dp1.10801[,keep_cols]
#dp1.10801.merge$geneticSampleID <- dp1.10801$geneticSampleID
merged <- merge(dp1.10801.merge, merged, all.x=T)


# Subset to 2014. we're not going to subset to peakGreenness for these covariates.
#merged <- merged[merged$sampleTiming == 'peakGreenness',]
merged$year <- substring(merged$dateID,1,4)
merged <- merged[merged$year == '2014',]
merged$siteID <- as.character(merged$siteID)
sites <- unique(merged$siteID)
# year_grab <- list()
# for(i in 1:length(sites)){
#   year_grab[[i]] <- min(merged[merged$siteID == sites[i], ]$year)
# }
# site_years <- data.frame(sites, unlist(year_grab))
# colnames(site_years) <- c('siteID','year')
# to_keep <- list()
# for(i in 1:nrow(site_years)){
#   to_keep[[i]] <- merged[merged$siteID == site_years$siteID[i] & merged$year == site_years$year[i],]
# }
# merged <- do.call(rbind, to_keep)
# # get rid of duped geneticSampleIDs.
merged <- merged[!(duplicated(merged)),]

# convert pH units from KCl to CaCl using the one-step converting function.
# this takes a dataset with both KCl and CaCl measurements, and uses their linear relationship for conversion
pH_conversion <- as.matrix(pH_conversion[,c("phca_val_1", "phkc_val_1")])
pH_KCl <- 
  convert_K(merged$soilInCaClpH, pH_conversion)
merged$pH_KCl <- pH_KCl$mean
merged$pH_sd <- pH_KCl$sd

# keep pH_water too, for bacterial forecasts
merged$pH_water <- merged$soilInWaterpH

#finalize columns for core.level.----
core.level <- merged[,c('sampleID','geneticSampleID','deprecatedVialID', 'dnaSampleID','siteID','plotID','dateID','collectDate','horizon','adjElevation','soilMoisture','sampleTiming','pH_water','pH_sd','pH_KCl','organicCPercent','CNratio')]
colnames(core.level)[(ncol(core.level) - 2) : ncol(core.level)] <- c('pH','pC','cn')
core.level$pC_sd <- pC_uncertainty_neon(core.level$pC)
core.level$cn_sd <- cn_uncertainty_neon(core.level$cn)

#get higher level observations and uncertainties (plot, site and global level.)----
#pC
pC.ag <- hierarch_core.means_JAGS(core.level$pC,core_plot = core.level$plotID)
pC.plot <- pC.ag$plot.table[,c('plotID','Mean','SD')]
pC.site <- pC.ag$site.table[,c('siteID','Mean','SD')]
pC.glob <- pC.ag$glob.table[,c('Mean','SD')]
colnames(pC.plot)[2:3] <- c('pC','pC_sd')
colnames(pC.site)[2:3] <- c('pC','pC_sd')
#C:N
cn.ag <- hierarch_core.means_JAGS(core.level$cn,core_plot = core.level$plotID)
cn.plot <- cn.ag$plot.table[,c('plotID','Mean','SD')]
cn.site <- cn.ag$site.table[,c('siteID','Mean','SD')]
cn.glob <- cn.ag$glob.table[,c('Mean','SD')]
colnames(cn.plot)[2:3] <- c('cn','cn_sd')
colnames(cn.site)[2:3] <- c('cn','cn_sd')
#pH - KCl (default)
pH.ag <- hierarch_core.means_JAGS(core.level$pH,core_plot = core.level$plotID)
pH.plot <- pH.ag$plot.table[,c('plotID','Mean','SD')]
pH.site <- pH.ag$site.table[,c('siteID','Mean','SD')]
pH.glob <- pH.ag$glob.table[,c('Mean','SD')]
colnames(pH.plot)[2:3] <- c('pH','pH_sd')
colnames(pH.site)[2:3] <- c('pH','pH_sd')

#pH - water (for bacteria)
pH_water.ag <- hierarch_core.means_JAGS(core.level$pH_water,core_plot = core.level$plotID)
pH_water.plot <- pH_water.ag$plot.table[,c('plotID','Mean','SD')]
pH_water.site <- pH_water.ag$site.table[,c('siteID','Mean','SD')]
pH_water.glob <- pH_water.ag$glob.table[,c('Mean','SD')]
colnames(pH_water.plot)[2:3] <- c('pH_water','pH_water_sd')
colnames(pH_water.site)[2:3] <- c('pH_water','pH_water_sd')

# now that plot- and site-level averages have used any available data from 2014,
# let's remove samples that are not from peakGreenness
core.level <- core.level[which(core.level$sampleTiming=="peakGreenness"),]

# let's also remove samples without a deprecatedVialID, which is how our samples are identified.
core.level <- core.level[!is.na(core.level$deprecatedVialID) & core.level$deprecatedVialID!="",]

#final output, core level.
core.out <- core.level

# removing "-GEN" from data.
core.out$geneticSampleID <- gsub('-GEN','',core.out$geneticSampleID)
dp1.10801$geneticSampleID <- gsub('-GEN','',dp1.10801$geneticSampleID)

#y obs need to be in core.out, and vice versa. 
# Essentially dropping any observation with no associated core-level covariate.
core.obs <- dp1.10801[dp1.10801$geneticSampleID %in% core.out$geneticSampleID,]
core.obs <- core.obs[!is.na(core.obs$deprecatedVialID) & core.obs$deprecatedVialID!="",]
core.obs <- core.obs[!(duplicated(core.obs$geneticSampleID)),] #drop duplicated, as before.
core.out <- core.out[core.out$geneticSampleID %in% core.obs$geneticSampleID,]
core.out <- core.out[!(duplicated(core.out$geneticSampleID)),] #drop duplicated, as before.

#plot level.
plot.out <- merge(pC.plot , cn.plot, all = T)
plot.out <- merge(plot.out, pH.plot, all = T)
plot.out <- merge(plot.out, pH_water.plot, all = T)
#Some plots in core-level not in plot-level.
to_add <- as.character(unique(core.out[!(core.out$plotID %in% plot.out$plotID),]$plotID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'plotID'
plot.out <- plyr::rbind.fill(plot.out,to_add)
plot.out$siteID <- substring(plot.out$plotID, 1, 4)
plot.out <- plot.out[order(plot.out$plotID),]
# remove plots in plot-level that aren't in core-level.
plot.out <- plot.out[plot.out$plotID %in% core.out$plotID,]

#site level.
site.out <- merge(pC.site, cn.site, all = T)
site.out <- merge(site.out, pH.site, all = T)
site.out <- merge(site.out, pH_water.site, all = T)
#Some sites in plot-level not in site-level.
to_add <- as.character(unique(plot.out[!(plot.out$siteID %in% site.out$siteID),]$siteID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'siteID'
site.out <- plyr::rbind.fill(site.out,to_add)
# remove sites in site-level that aren't in plot-level.
site.out <- site.out[site.out$siteID %in% plot.out$siteID,]

#global level
glob.out <- rbind(pC.glob, cn.glob, pH.glob, pH_water.glob)
pred <- c('pC','cn','pH','pH_water')
glob.out <- cbind(pred,glob.out)

#save output.----
saveRDS(core.obs,  core_obs_data.path)
saveRDS(core.out, core_core_data.path)
saveRDS(plot.out, core_plot_data.path)
saveRDS(site.out, core_site_data.path)
saveRDS(glob.out, core_glob_data.path)
