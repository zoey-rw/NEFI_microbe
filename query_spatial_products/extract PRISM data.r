#writing function to get mean annual temperature and precipitation from PRISM.
#Also extract that month's tmean, tmin, tmax and ppt.
#requires some packages
rm(list=ls())
library(sp)
library(raster)
source('NEFI_functions/prism_query.r')

#path to your data file
data.path <- '/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds'

#grab your data.
d <- readRDS(data.path)

#extract prism normals and month data
#prism path now set as default to directory in my data folder on pecan2.
climate_data <- prism_query(d)

#append this data to your dataframe, save it.
d <- cbind(d,climate_data)
saveRDS(d,data.path)