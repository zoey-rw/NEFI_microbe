#site_obs_aggregation.r
rm(list=ls())
source('paths.r')

site.out <- readRDS(site_level_data_16S.path)

#get global means.
to_ag <- site.out[,c('map','mat','NPP','n.dep','dry.dep')]
 Mean <- apply(to_ag, 2, mean)
   SD <- apply(to_ag, 2,   sd)
preds <- colnames(to_ag)
glob.out <- data.frame(preds,Mean, SD)

#Get forest and conifer 0-1 assignments.
siteID <- c("ORNL", "CPER", "WOOD", "TALL", "JERC", "OSBS", "RMNP", "HARV", "BART", "STER", "SCBI", "DSNY", "UNDE")
forest <- c(1,0,0,1,0,1,1,1,1,0,1,0,1)
conifer <- c(1,0,0,1,0,1,1,1,1,0,0,1,1)
moisture <- extract_soil_moist(site.out[,2], site.out[,3])
to_merge <- data.frame(siteID,forest,conifer,moisture)  
site.out <- merge(site.out, to_merge, all.x=T)

#save output.
saveRDS(site.out, site_site_16S.path)
saveRDS(glob.out, site_glob_16S.path)
