#Extract data from NEON.DP1.10078.001 - periodic measurements of soil C and N.

rm(list=ls())
source('paths.r')
library(neonUtilities)

# Specify output path
output.path <- dp1.10078.001_output.path

# Grab all sites that have soil core data
# Downloads about 1 MB, from 2013 to 2019, into the workspace. Takes a few minutes.
dat <- loadByProduct(dpID="DP1.10078.001", site="all", startdate="2013-01", enddate="2019-05",
                     package="expanded", check.size=F)

# Rename dataframes
output <- dat$sls_soilChemistry

# Add a couple columns for later
output$dateID <- substr(output$collectDate,1,7)
output$site_date <- paste0(output$siteID,'-', output$dateID)
output$site_date_plot <- paste0(output$siteID,'-',output$dateID,'-',output$plotID)


# Save DP.10078 output
saveRDS(output, output.path)
