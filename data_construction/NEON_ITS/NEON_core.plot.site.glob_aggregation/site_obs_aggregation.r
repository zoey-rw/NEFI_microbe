#site_obs_aggregation.r
rm(list=ls())
source('paths.r')

site.out <- readRDS(site_level_data.path)

#get global means.
to_ag <- site.out[,c('map','mat','NPP','n.dep','dry.dep')]
 Mean <- apply(to_ag, 2, mean)
   SD <- apply(to_ag, 2,   sd)
preds <- colnames(to_ag)
glob.out <- data.frame(preds,Mean, SD)

#save output.
saveRDS(site.out, site_site.path)
saveRDS(glob.out, site_glob.path)