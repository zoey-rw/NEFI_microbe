#Formatting_NEON_microbial - Colin Averill, October 2, 2017
#This script uses the Mapping files and OTU tables provided by Lee Stanish for NEON microbial data.
#It first subsets and orders mapping and OTU tables so that each contain the same observations.
#It then grabs site-level latitude and longitude for each site through the nneo package, and appends this to the mapping file.
#It then pulls apart the samples codes of the mapping file to get site and date information.
#Finally, this script saves these formatted OTU and mapping files for downstream querying of other NEON products.

#clear R environment, load packages
rm(list=ls())
library(data.table)
library(nneo)

#set output paths
output.dir <- "/fs/data3/caverill/NEFI_microbial/map_otu/" #output directory
output.map <- "16S_map_clean.rds" 
output.otu <- "16S_otu_clean.rds"
output.map <- paste0(output.dir,output.map)
output.otu <- paste0(output.dir,output.otu)

#subtring function to pull from the right edge of a character entry.
substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

#load raw mapping file and OTU tables as received from Lee Stanish.
map <- (read.csv('/fs/data3/caverill/NEFI_microbial/map_otu/16S_combined_mapping_qced_Summer_2017.csv'))
otu <- (read.table('/fs/data3/caverill/NEFI_microbial/map_otu/otu_table_16S.txt', header=T))


#First, make sure all map observations are in the OTU table, and vice versa. Discard those that are not.
map <- map[map$Description %in% colnames(otu),]
otu <- otu[,colnames(otu) %in% map$Description]

#sort map and otu table so they match.
map <- map[order(map$Description),]
otu <- otu[,order(colnames(otu))]

#Convert genetic sample IDs to uppercase so they match the NEON data products.
map$geneticSampleID <- toupper(map$geneticSampleID)

#Grab site IDs for each observation.
map$site <- substr(map$geneticSampleID,1,4)

#Get site-level latitude and longitude for each observation.
#grab unique sites
sites <- unique(map$site)
#make it a dataframe.
sites <- data.frame(sites)
sites$latitude  <- NA
sites$longitude <- NA

#loop over sites to extract latitude/longitude from nneo package.
for(i in 1:nrow(sites)){
  sites$latitude[i] <- nneo_location(sites$sites[i])$locationDecimalLatitude
  sites$longitude[i] <- nneo_location(sites$sites[i])$locationDecimalLongitude
}
#merge lat-long into mapping file
map <- merge(map,sites,by.x='site',by.y='sites',all.x=T)

#pull apart the genetic Sample ID to get the date. 
map$geneticSampleID <- as.character(map$geneticSampleID)
map$date <- substrRight(map$geneticSampleID,12)
map$date <- substr(map$date,1,8)
map$year  <- substr(map$date,1,4)
map$month <- substr(map$date,5,6)
map$day   <- substr(map$date,7,8)
map$date <- paste0(map$year,'-',map$month,'-',map$day)
map$year_month <- paste0(map$year,'-',map$month)
#Give it a continuous date value, with 0 = Jan 1, 2013
date.lookup <- format(seq(as.Date("2013-01-01"), as.Date("2016-12-31"), by = "1 day"))
map$epoch_date <- match(map$date, date.lookup)

#Save output
saveRDS(map,output.map)
saveRDS(otu,output.otu)