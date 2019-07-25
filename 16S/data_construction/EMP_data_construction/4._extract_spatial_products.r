#Assign spatial products.
#Climate from worldclim2. N-dep from NADP.
#These are based on functions written by Colin. 
#This script will only run on pecan2 as this is where these spatial products are actually stored (they are big).
rm(list = ls())
library(data.table)
source('paths.r')
source('NEFI_functions/extract_ndep.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')
source('NEFI_functions/extract_pH.r')
source('NEFI_functions/extract_npp.r')


#load 16S meta data file with lat/long
d <- readRDS(emp_map_clean.path)

#make sure things you are going to extract aren't already in dataframe. If they are, remove them so we dont double name things.
to_remove <- c('mat','map','mat_sd','map_sd','mat_CV','map_CV','mdr','aridity','n.dep','dry.dep','wet.dep','ph_estim','ndep.glob')
d <- d[,!(colnames(d) %in% to_remove)]

#extract worldclim.
clim <- worldclim2_grab(d$latitude, d$longitude, d$elevation_m, worldclim2_folder = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/WorldClim2/")
clim$aridity <- arid_extract(d$latitude,d$longitude, folder = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/Global_Aridity/")

#extract N deposition. Half of these are NA because our Ndep product only covers the United States.
#ignore the warnings.
ndep <- extract_ndep(d$longitude, d$latitude, folder = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/CASTNET_Ndep/")

#extract pH
ph_estim <- extract_pH(d$latitude, d$longitude, folder = "/projectnb/talbot-lab-data/spatial_raster_data/SoilGrids_uncertainty/")

NPP     <- extract_npp(d$latitude, d$longitude, path="/projectnb/talbot-lab-data/NEFI_data/covariate_data/NPP/w001001.adf")

# from Ackerman et al.: https://conservancy.umn.edu/handle/11299/197613
ndep.glob <- read.csv("/projectnb/talbot-lab-data/NEFI_data/covariate_data/inorganic_N_deposition.csv")
ndep.tot.2014 <- n.dep[,c("latitude", "longitude", "pixel_area_km2","tot_2014")]
ndep.global <- list()
for (n in 1:nrow(d)){
  lat <- round(d$latitude)[n]
  lon <- round(d$longitude)[n]
  sorted_lat <- sort(ndep.tot.2014$latitude)
  newlat <- sorted_lat[MALDIquant::match.closest(lat, sorted_lat)]
  sorted_lon <- sort(ndep.tot.2014$longitude)
  newlon <- sorted_lon[MALDIquant::match.closest(lon, sorted_lon)]
  
  ndep.point <- n.dep.tot.2014[which(ndep.tot.2014$latitude==newlat & ndep.tot.2014$longitude==newlon),]
  ndep.global[n] <- ndep.point$tot_2014/100 #convert per-km to per-hectare
}  

#put together all spatial products
spatial <- cbind(ndep, ph_estim,clim, NPP, ndep.glob = unlist(ndep.global))

#remove any columns from previous spatial extractions, update with new extraction.
d <- as.data.frame(d)
d <- cbind(d,spatial)

#save
saveRDS(d,emp_map_clean.path)
