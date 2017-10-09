#2._extracting_soil_physical_properties.r Colin Averill, October 3, 2017
#This script extracts physical properties for all site by month combinations from a mapping file.
#This will be generalized to a function that takes as input a vector of sites and year-months.

#clear R environment, load packages.
rm(list=ls())
library(nneo)
library(data.table)
#source the function that queries soil physical properties.
source('NEFI_functions/soil_physical_query.r')

#grab mapping file of interest.
map <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')

#use function that queries the neon soil physical properties data product.
out.physical <- soil_physical_query(map)

#merge together map and physical soil attributes by genetic sample ID, cause thats what Colin cares about.
#all.x=T tells merge to enter NAs if those data products don't have data.
test <- merge(map, out.physical, 
              all.x = T, 
              by='geneticSampleID')
