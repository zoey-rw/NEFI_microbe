#Extract data from NEON.DP1.10098.001 - NEON forest inventory data.
#clear environment, source output paths.
rm(list=ls())
source('paths.r')
library(neonUtilities)

# Specify output path
output.path <- dp1.10098.001_output.path

# Grab all sites that have soil core data
# Downloads about 70 MB, from 2013 to 2019, into the workspace. Takes ~10 minutes.
dat <- loadByProduct(dpID="DP1.10098.001", site="all", startdate="2013-01", enddate="2019-05",
                     package="expanded", check.size=F)

# Rename dataframes
diam <- dat$vst_apparentindividual # individual diameters 
spp <- dat$vst_mappingandtagging # species IDs

# Merge species info into individual tree info
spp.merge <- spp[,!(colnames(spp) %in% colnames(diam))]
spp.merge$individualID <- spp$individualID
output <- merge(diam,spp.merge, all = T)

# Add a couple columns for later
output$dateID <- substr(output$date,1,7)
output$site_date <- paste(output$siteID, output$dateID)

# Save DP1.10098 output
saveRDS(output, output.path)
