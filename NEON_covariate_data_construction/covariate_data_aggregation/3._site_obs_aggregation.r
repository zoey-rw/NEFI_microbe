#site_obs_aggregation.r
rm(list=ls())
source('paths.r')

# load site-level data and core-level data
site.out <- readRDS(site_level_covariates.path)
core_obs <- readRDS(core_obs_data.path)

# subset to the sites that we have observations for
site.out <- site.out[site.out$siteID %in% core_obs$siteID,]

# convert MAP to from mm to m.
site.out$map <- site.out$map/1000
site.out$map_sd <- site.out$map_sd/1000

#get global means.
to_ag <- site.out[,c('map','mat','NPP','n.dep','dry.dep','ndep.glob')]
Mean <- apply(to_ag, 2, mean, na.rm=TRUE)
SD <- apply(to_ag, 2,   sd, na.rm=TRUE)
preds <- colnames(to_ag)
glob.out <- data.frame(preds,Mean, SD)

#Get forest and conifer 0-1 assignments.
siteID <- c("ORNL", "CPER", "WOOD", "TALL", "JERC", "OSBS", "RMNP", "HARV", "BART", "STER", "SCBI", "DSNY", "UNDE")
forest <- c(1,0,0,1,0,1,1,1,1,0,1,0,1)
conifer <- c(1,0,0,1,0,1,1,1,1,0,0,1,1)
#moisture <- extract_soil_moist(site.out[,2], site.out[,3])
to_merge <- data.frame(siteID,forest,conifer)  
site.out <- merge(site.out, to_merge, all.x=T)

#save output.
saveRDS(site.out, site_site_data.path)
saveRDS(glob.out, site_glob_data.path)
