#Plot level aggregation.
#Getting plot, site and global means and sd's for all observations made at the plot level.
#This is currently total basal area and %ECM trees.
#clear environment, source paths.
rm(list=ls())
source('paths.r')
source('NEFI_functions/hierarch_plot.means_JAGS.r')

#load data
plot_data <- readRDS(dp1.10098.00_plot.level.path)

#y-values to constrain some things.
core_obs <- readRDS(core_obs.path)

#get plot level uncertainty on relEM
#if relEM <= 0, logit transform to -10.
#if relEM => 1, logit transform to  10.
relEM.sim <- list()
for(i in 1:1000){
  all <- rnorm(nrow(plot_data),plot_data$basal_live, plot_data$basal_live_sd)
   EM <- rnorm(nrow(plot_data),plot_data$basal_ECM , plot_data$basal_ECM_sd)
   relEM <- EM / all
   for(k in 1:length(relEM)){
     if(relEM[k] < 0.0001){relEM[k] <- 0.0001}
     if(relEM[k] > 0.9999){relEM[k] <- 0.9999}
   }
   b.relEM <- boot::logit(relEM)
   relEM.sim[[i]] <- b.relEM
}
relEM.sim <- do.call(cbind,relEM.sim)
b.relEM <- rowMeans(relEM.sim)
b.relEM_sd <- matrixStats::rowSds(relEM.sim)

#plot level output.
plot.out <- plot_data
plot.out$b.relEM <- b.relEM
plot.out$b.relEM_sd <- b.relEM_sd

#aggregate up to site and glboal level
basal_live_ag <- hierarch_plot.means_JAGS(plot_data$basal_live, plot_site = plot_data$siteID)
basal_dead_ag <- hierarch_plot.means_JAGS(plot_data$basal_dead, plot_site = plot_data$siteID)
basal_ECM_ag  <- hierarch_plot.means_JAGS(plot_data$basal_ECM , plot_site = plot_data$siteID)
basal_bECM_ag <- hierarch_plot.means_JAGS(plot.out$b.relEM   , plot_site = plot.out$siteID)

#add plots to plot-level that are in y-observations but not present in this data product yet.
to_add <- as.character(unique(core_obs[!(core_obs$plotID %in% plot.out$plotID),]$plotID))
to_add.site <- substring(to_add,1,4)
to_add <- data.frame(to_add.site, to_add)
colnames(to_add) <- c('siteID','plotID')
plot.out <- plyr::rbind.fill(plot.out, to_add)

#site and global level output.
basal_live.site <- basal_live_ag$site.table[,c('siteID','Mean','SD')]
basal_live.glob <- basal_live_ag$glob.table[,c('Mean','SD')]
colnames(basal_live.site)[2:3] <- c('basal_live','basal_live_sd')
basal_dead.site <- basal_dead_ag$site.table[,c('siteID','Mean','SD')]
basal_dead.glob <- basal_dead_ag$glob.table[,c('Mean','SD')]
colnames(basal_dead.site)[2:3] <- c('basal_dead','basal_dead_sd')
basal_ECM.site  <- basal_ECM_ag$site.table [,c('siteID','Mean','SD')]
basal_ECM.glob  <- basal_ECM_ag$glob.table [,c('Mean','SD')]
colnames(basal_ECM.site)[2:3] <- c('basal_ECM','basal_ECM_sd')
basal_bECM.site <- basal_bECM_ag$site.table[,c('siteID','Mean','SD')]
basal_bECM.glob <- basal_bECM_ag$glob.table[,c('Mean','SD')]
colnames(basal_bECM.site)[2:3] <- c('b.relEM','b.relEM_sd')

#site level merge.
site.out <- merge(basal_live.site, basal_dead.site, all = T)
site.out <- merge(site.out, basal_ECM.site, all = T)
site.out <- merge(site.out, basal_bECM.site, all = T)
#add in sites in y data not present here.
to_add <- as.character(unique(core_obs[!(core_obs$siteID %in% site.out$siteID),]$siteID))
to_add <- data.frame(to_add)
colnames(to_add) <- 'siteID'
site.out <- plyr::rbind.fill(site.out,to_add)

#global level merge.
glob.out <- rbind(basal_live.glob,basal_dead.glob,basal_ECM.glob, basal_bECM.glob)
pred <- c('basal_live','basal_dead','basal_ECM','b.relEM')
glob.out <- cbind(pred,glob.out)

#save output.----
saveRDS(plot.out, plot_plot.path)
saveRDS(site.out, plot_site.path)
saveRDS(glob.out, plot_glob.path)

