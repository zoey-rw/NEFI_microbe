#Get soil physical properties data from DP1.10086.00 for 16S and ITS data

rm(list=ls())
source('paths.r')
library(geoNEON)
library(neonUtilities)

# Specify output path
output.path <- dp1.10086.001_output.path

# Grab all sites that have soil core data
# Downloads about 15 MB, from 2013 to 2019, into the workspace. Takes a few minutes.
dat <- loadByProduct(dpID="DP1.10086.001", site="all", startdate="2013-01", enddate="2019-05",
                     package="expanded", check.size=F)

# Rename dataframes
core.data <- dat$sls_soilCoreCollection
moisture.data <- dat$sls_soilMoisture
pH.data <- dat$sls_soilpH

# Merge moisture into soil cores
moist.merge <- moisture.data[,!(colnames(moisture.data) %in% colnames(core.data))]
moist.merge$sampleID <- moisture.data$sampleID
core.data <- merge(core.data,moist.merge, all = T)

# Merge pH into soil cores
pH.merge <- pH.data[,!(colnames(pH.data) %in% colnames(core.data))]
pH.merge$sampleID <- pH.data$sampleID
output <- merge(core.data,pH.merge, all = T)

# Add a couple columns for later
output$dateID <- substr(output$collectDate,1,7)
output$site_date <- paste(output$siteID, output$dateID)

# Add precise geolocation data
output <- geoNEON::def.calc.geo.os(output, 'sls_soilCoreCollection')

# Save DP1.10086 output
saveRDS(output, output.path)