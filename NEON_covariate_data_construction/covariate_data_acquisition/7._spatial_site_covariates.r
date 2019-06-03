# Get site level climate and N-deposition data with associated uncertainties where applicable.
# clear environment, load paths and custom functions.
# ignore the warnings on the N-deposition extraction.
rm(list = ls())
source('paths.r')
source('NEFI_functions/arid_extract.r')
source('NEFI_functions/extract_ndep.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/extract_npp.r')

# Specify output path
output.path <- site_level_covariates.path

# Get site level lat-longitude.
location <- readRDS(dp1.10086.001_output.path)
sites <- unique(location$siteID)

core_obs <- readRDS(core_obs_data.path)
# subset to the sites that we have observations for
sites <- sites[sites %in% core_obs$siteID]
location <- location[location$siteID %in% sites,]

# Within eah site get mean latitude, longitude and elevation.
lat <- aggregate(adjDecimalLatitude  ~ siteID, data = location, FUN = mean)
lon <- aggregate(adjDecimalLongitude ~ siteID, data = location, FUN = mean)
ele <- aggregate(adjElevation        ~ siteID, data = location, FUN = mean)

# Extract climate and uncertainty from worldclim2
climate <- worldclim2_grab(latitude = lat[,2], longitude = lon[,2], elev = ele[,2])
climate$aridity <- arid_extract(lat[,2],lon[,2])
# Extract wet, dry and total nitrogen deposition.
ndep    <- extract_ndep(longitude = lon[,2],latitude = lat[,2])
# Extract NPP
NPP     <- extract_npp(latitude = lat[,2], longitude = lon[,2])

# Wrap up output and save.
output <- cbind(lat,lon[,2],ele[,2],climate,ndep,NPP)
colnames(output)[1:4] <- c('siteID','latitude','longitude','elevation')
saveRDS(output,output.path)