##### preparing Bahram prior metadata. #####
# later combined with abundance data for the prior model fits.

##### clear environment, source packages, functions and paths #####
rm(list=ls())

# devtools::install_version("lubridate", "1.6.0") # something about Red Hat Enterprise Linux 6; currently need to use older version
library(lubridate)  
library(data.table)
library(rgdal)
library(raster)

#library(SRAdb)
source('paths.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')

# set output path
output.path <- bahram_metadata.path 


##### load files ####

# SRA RunInfo table from: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP021922
SRA <- read.csv(SRA.path)

# metadata file that was sent to Colin
metadata_raw <- read.csv(metadata_bahram_raw.path)

# load tedersoo mapping file
map <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))
map <- data.table(map)

# load times - sent to Colin separately by Leho Tedersoo
time <- read.csv(ted_sampling_dates.path, header = TRUE, row.names=1, check.names = FALSE)


##### Format metadata file #####

# combine SRA sample names with Bahram sample names.
sample_names <- SRA[,c(7,9,11)]
names(sample_names)[names(sample_names) == 'environment_biome'] <- 'Biome'
metadata <- merge(sample_names, metadata_raw, by.x = "Sample_Name", by.y = "X")
metadata <- data.table(metadata)

# subset to temperate latitudes
metadata <- metadata[Lat < 66.5 & Lat > 23.5,]

# merge in site, moisture, and relEM from Tedersoo file
metadata <- metadata[,-c("Moisture")]
metadata <- merge(x = metadata, 
                  y = map[ , c("tedersoo.code", "Site", "Moisture", 
                               "Relative.basal.area.of.EcM.trees.....of.total.basal.area.of.all.AM.and.EcM.tees.taken.together.")], 
                  by.x = "Sample_Name", by.y="tedersoo.code", all.x=TRUE)

# format the time data frame (get rid of an empty column, etc.)
colnames(time)[1] <- 'human.date'
time[2] <- NULL
#convert human readable date to days since epoch
time$epoch.date <- round(
  as.numeric(
    as.POSIXlt(time$human.date, 
               format = "%m/%d/%y", origin = "01/01/1970"))/86400)
# get day of year (doy) as well.
time$doy  <- lubridate::yday(as.Date(time$human.date,format='%m/%d/%Y'))
time$year <- lubridate::year(as.Date(time$human.date,format='%m/%d/%Y'))
time$epoch.date <- (time$year - 9)*365 + time$doy

# drop samples in time table not present in mapping file.
time <- time[row.names(time) %in% metadata$Sample_Name,]

# push times into the mapping file
time$Sample_Name <- rownames(time)
metadata <- merge(metadata,time, by = 'Sample_Name', all.x=T)

# assign forest variables
# shouldn't this also include  'Moist tropical forests', 'Tropical montane forests' and 'Southern_temperate_forests'?
metadata$forest <-ifelse(metadata$Biome %in% c('Temperate_coniferous_forests',
                                               'Temperate_deciduous_forests','Dry_tropical_forests','Boreal_forests'),1,0)
metadata$conifer <- ifelse(metadata$Biome %in% c('Temperate_coniferous_forests'),1,0)



#get worldclim2 climate variables and aridity index, merge with metadata
climate <- worldclim2_grab(latitude = metadata$Lat,longitude=metadata$Lon)
climate$aridity <- arid_extract(metadata$Lat, metadata$Lon)
metadata <- cbind(metadata, climate)

# grab columns from map actually of interest.
metadata_subset <- metadata[,.(Sample_Name,Run,Lon,Lat,PH,Moisture,N,C,C.N,human.date,doy,epoch.date,NPP,forest,conifer)]

# rename columns
setnames(metadata,c('Sample_Name','Moisture','N' ,'C' ,'C.N',
                       'Relative.basal.area.of.EcM.trees.....of.total.basal.area.of.all.AM.and.EcM.tees.taken.together.'),
         c('Mapping.ID','moisture','pN','pC','cn','relEM'))

# save output.----
saveRDS(metadata, output.path)
