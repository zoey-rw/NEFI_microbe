#Querying N-deposition data from NADP and CASTNET products from the EPA.
#clear environment, load packages, source functions.
rm(list=ls())
library(rgdal)
library(sp)
library(raster)
source('NEFI_functions/extract_ndep.r')

#load data file with lat/long
d <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')

out <- extract_ndep(d$longitude,d$latitude)

d <- cbind(d,out)
