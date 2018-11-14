##### preparing Bahram prior, using taxonomic table and metadata. #####

# outputs 20 most cosmopolitan genera
# functional groups are not yet assigned

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

# number of genera to keep (top 10 or 20 most abundant)
n.gen <- 20

##### load files ####
dir <- data.dir 
# SRA RunInfo table from: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP021922
SRA <- read.csv(paste0(dir,"16S/scc_gen/SraRunTable.csv"))

# metadata file that was sent to Colin
metadata_raw <- read.csv(paste0(dir,"16S/scc_gen/metadata_bahram.csv"))

# load tedersoo mapping file
map <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))
map <- data.table(map)

# load SV table as otu file
otu <- readRDS(bahram_dada2_SV_table.path)

# load taxonomy
tax <- readRDS(bahram_dada2_tax_table.path)

# load times - sent to Colin separately by Leho Tedersoo
time <- read.csv(ted_sampling_dates.path, header = TRUE, row.names=1, check.names = FALSE)


##### Format metadata file #####

# combine SRA sample names with Bahram sample names.
sample_names <- SRA[,c(7,9,11)]
names(sample_names)[names(sample_names) == 'environment_biome'] <- 'Biome'
metadata <- merge(sample_names, metadata_raw, by.x = "Sample_Name", by.y = "X")
metadata <- data.table(metadata)

# subset to northern temperate latitudes
metadata <- metadata[Lat < 66.5 & Lat > 23.5,]

# merge in site and moisture from Tedersoo file
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

##### taxonomic and functional assignment #####

# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# subset otu table and tax table to only include observations in map file
metadata$Run <- as.character(metadata$Run)
otu <- otu[rownames(otu) %in% metadata$Run,]
metadata <- metadata[metadata$Run %in% rownames(otu),]
# order OTU table to match the mapping file
otu <- otu[order(rownames(otu), metadata$Run),]

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'Bacteria'|tax$kingdom == 'Archaea',] # only removes ~500 counts
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]



#determine cosmopolitan genera.----
#condition 1: present in greater than 50% of samples.
#condition 2: top 10% sequence abundance (sensu Delgado-Barquez et al. 2018, Science)
genera <- unique(tax$genus)
genera <- as.character(genera)
test <- data.table(cbind(tax, t(otu)))
seq.out <- list()
cosmo.out <- list()
for(i in 1:length(genera)){
  z <- test[genus == genera[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  cosmo <- length(out[out > 0]) / length(out)
  seq.out[[i]] <- out
  cosmo.out[[i]] <- cosmo
}
seq.out <- do.call('rbind',seq.out)
cosmo.out <- do.call('rbind',cosmo.out)
j <- data.table(cbind(genera,cosmo.out))
colnames(j)[2] <- 'cosmo'
j$cosmo <- as.numeric(j$cosmo)
j <- j[order(-cosmo),]
counts <- rowSums(seq.out)
k <- data.table(cbind(genera,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(genera %in% c('unidentified'))] #remove the genus "unidentified".
k <- k[genera!=""&!is.na(genera),] #remove NA and empty genera
j <- j[!(genera %in% c('unidentified'))] #remove the genus "unidentified".
j <- j[genera!=""&!is.na(genera),] #remove NA and empty genera
head(cbind(j,k), 20)

#7 genera present in > 50% of samples are also in the top 10%.
#cosmo_genera <- j[cosmo >= 0.50,]$genera

# 20 genera present in > 94% of samples are also in the top 10%.
# maybe? where is the 10% calculated?
cosmo_genera <- j[cosmo >= 0.94,]$genera

#Get seq abundances of cosmo genera.----
gen.list <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(cosmo_genera)){
  z <- k[genus == cosmo_genera[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  gen.list[[i]] <- out
}
gen.list <- data.frame(t(do.call('rbind',gen.list)))
colnames(gen.list) <- cosmo_genera
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(gen.list)
gen.list <- cbind(other,gen.list)
gen.list <- list(gen.list,seq_total)
names(gen.list) <- c('abundances','seq_total')
gen.list$rel.abundances <- gen.list$abundances / gen.list$seq_total


#save output.----
saveRDS(gen.list, cosmo_output_16S.path)










# don't think we need all of this anymore - remove genera in below script, just needs to be metadata
# 
# # normalize the otu table 
# otu <- t(otu)
# pro.function <- function(otu){
#   for(i in 1:ncol(otu)){
#     otu[,i] <- otu[,i] / sum(otu[,i])
#   }
#   return(otu)
# }
# otu <- pro.function(otu)
# 
# # make sure column sums are 1
# colSums(otu)
# 
# 
# 
# #### get most abundant genera in dataset.####
# genera <- unique(tax$genus)
# test <- data.table(cbind(tax, otu))
# seq.out <- list()
# for(i in 1:length(genera)){
#   z <- test[genus == genera[i],]
#   start <- ncol(tax) + 1
#   out <- colSums(z[,start:ncol(z)])
#   seq.out[[i]] <- out
# }
# 
# # Count genus level abundance, grab some number of most abundant genera
# seq.out <- do.call('rbind',seq.out)
# counts <- rowSums(seq.out)
# genera <- as.character(genera)
# k <- data.table(cbind(genera,counts))
# k$counts <- as.numeric(as.character(k$counts))
# k <- k[order(-counts),]
# k <- k[genera!=""&!is.na(genera),] #remove NA and empty genera
# #grab genera of interest.
# of_interest <- k$genera[1:n.gen]
# 
# # Get relative abundances of the most abundant genera
# gen.list <- list()
# for(i in 1:length(of_interest)){
#   z <- data.table(cbind(tax,otu))
#   z <- z[genus %in% of_interest[i],]
#   start <- ncol(tax) + 1
#   out <- colSums(z[,start:ncol(z)])
#   gen.list[[i]] <- out
# }
# gen.list <- data.frame(t(do.call('rbind',gen.list)))
# colnames(gen.list) <- of_interest
# gen.list$Mapping.ID <- rownames(gen.list)
# 
# 
# #### get most cosmopolitan genera in dataset ####
# 
# genera <- unique(tax$genus)
# test <- data.table(cbind(tax, otu))
# seq.out <- list()
# for(i in 1:length(genera)){
#   z <- test[genus == genera[i],]
#   start <- ncol(tax) + 1
#   out <- colSums(z[,start:ncol(z)])
#   seq.out[[i]] <- out
# }
# 
# # Count genus level abundance
# seq.out <- do.call('rbind',seq.out)
# counts <- rowSums(seq.out)
# genera <- as.character(genera)
# k <- data.table(cbind(genera,counts))
# k$counts <- as.numeric(as.character(k$counts))
# k <- k[order(-counts),]
# k <- k[genera!=""&!is.na(genera),] #remove NA and empty genera
# #grab genera of interest (all).
# of_interest <- k$genera 
# 
# # Get relative abundances of all genera
# gen.list <- list()
# for(i in 1:length(of_interest)){
#   z <- data.table(cbind(tax,otu))
#   z <- z[genus %in% of_interest[i],]
#   start <- ncol(tax) + 1
#   out <- colSums(z[,start:ncol(z)])
#   gen.list[[i]] <- out
# }
# gen.list <- data.frame(t(do.call('rbind',gen.list)))
# colnames(gen.list) <- of_interest
# gen.list$Mapping.ID <- rownames(gen.list)
# 
# all.genera.abundances <- gen.list
# sums <- colSums(all.genera.abundances == 0)
# cosmo.gen <- names(sums[sums<2])
# cosmo.gen <- cosmo.gen[cosmo.gen != "Mapping.ID"]
# 
# 
# # Get relative abundances of the most cosmopolitan genera
# of_interest <- cosmo.gen
# gen.list <- list()
# for(i in 1:length(of_interest)){
#   z <- data.table(cbind(tax,otu))
#   z <- z[genus %in% of_interest[i],]
#   start <- ncol(tax) + 1
#   out <- colSums(z[,start:ncol(z)])
#   gen.list[[i]] <- out
# }
# gen.list <- data.frame(t(do.call('rbind',gen.list)))
# colnames(gen.list) <- of_interest
# gen.list$Mapping.ID <- rownames(gen.list)
# cosmo.gen.abundances <- gen.list
# 
# 
# 
# 
# 
# 
# 
# # Merge together data aggregated groups you want for analysis.
# # abundances <- merge(fun.list,gen.list)
# abundances <- cosmo.gen.abundances
# 
# #get worldclim2 climate variables and aridity index, merge with metadata
# climate <- worldclim2_grab(latitude = metadata$Lat, longitude = metadata$Lon)
# climate$aridity <- arid_extract(metadata$Lat, metadata$Lon)
# metadata <- cbind(metadata, climate)
# 
# # final merging of abundances and metadata
# # grab columns from map actually of interest.
# metadata_subset <- metadata[,.(Sample_Name,Run,Lon,Lat,PH,Moisture,N,C,C.N,human.date,doy,epoch.date,NPP,forest,conifer)]
# metadata <- merge(metadata,abundances, by.x = 'Run',by.y = 'Mapping.ID')
# 
# # rename columns
# setnames(metadata,c('Sample_Name','Moisture','N' ,'C' ,'C.N', 
#                        'Relative.basal.area.of.EcM.trees.....of.total.basal.area.of.all.AM.and.EcM.tees.taken.together.'), 
#          c('Mapping.ID','moisture','pN','pC','cn','relEM'))
# 
# # save output.----
# saveRDS(metadata, bahram_prior.path)
