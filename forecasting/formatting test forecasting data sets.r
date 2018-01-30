#generate some forecasting data
rm(list=ls())
library(betareg)
library(raster)
library(sp)
library(nneo)
library(data.table)
library(lubridate)
source('NEFI_functions/prism_query.r')
source('NEFI_functions/extract_ndep.r')
source('NEFI_functions/soil_physical_query.r')
source('NEFI_functions/neon_query.r')

#set output paths for core level data
o.dir <- '/fs/data3/caverill/NEFI_microbial/prior_data/'
 tal.core.path <- paste0(o.dir,'tal_core_data.rds')
 tal.site.path <- paste0(o.dir,'tal_site_data.rds')
neon.core.path <- paste0(o.dir,'neon_core_data.rds')
neon.site.path <- paste0(o.dir,'neon_site_data.rds')

#load Talbot data for prior and neon fungal data for forecast.
tal <- data.table(readRDS('/fs/data3/caverill/Microbial_Space_Time_data/talbot_2014.data/tal_map_filtered.rds'))
neon.map <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/ITS_map_clean.rds')
neon.otu <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/ITS_otu_clean.rds')
neon.tax <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/ITS_tax_clean.rds')

#subset talbot data
tal.prior <- tal[,.(Site,Plot,Horizon,longitude.dd,latitude.dd,Perc.N,Perc.C,Perc.Soil.Moisture,pH,Bulk.Density.,Date.Sampled,em.seqs,SAP,AMF,White_rot,total.seqs,MAT,MAP)]
tal.prior$relEM <- tal.prior$em.seqs / tal.prior$total.seqs
tal.prior$relAM <- tal.prior$AMF / tal.prior$total.seqs
tal.prior$relSAP <- tal.prior$SAP / tal.prior$total.seqs
tal.prior$relWR <- tal.prior$White_rot / tal.prior$total.seqs
tal.prior$cn <- tal.prior$Perc.C / tal.prior$Perc.N

#get date together
tal.prior$doy <- as.Date(tal.prior$Date.Sampled,format='%m/%d/%Y')
tal.prior$year <- lubridate::year(tal.prior$doy)
tal.prior$year <- paste0('20',tal.prior$year)
tal.prior$month <- lubridate::month(tal.prior$doy)
tal.prior$month <- formatC(tal.prior$month, width = 2, flag = "0") # add leadingzero to month.
tal.prior$day  <- lubridate::day(tal.prior$doy)
tal.prior$day  <- formatC(tal.prior$day, width = 2, flag = "0") # add leadingzero to day
tal.prior$doy  <- lubridate::yday(tal.prior$doy)
tal.prior$year_month <- paste(tal.prior$year, tal.prior$month, sep = '_')
tal.prior$date <- paste(tal.prior$year, tal.prior$month, tal.prior$day, sep = "-")
#Give it a continuous date value, with 0 = Jan 1, 2010
date.lookup <- format(seq(as.Date("2010-01-01"), as.Date("2016-12-31"), by = "1 day"))
tal.prior$epoch_date <- match(tal.prior$date, date.lookup)

#change talbot prior column names to match neon column names.
setnames(tal.prior,'longitude.dd','longitude')
setnames(tal.prior,'latitude.dd','latitude')
setnames(tal.prior,'Horizon','horizon')
setnames(tal.prior,'MAP','map30')
setnames(tal.prior,'MAT','mat30')
setnames(tal.prior,'Site','site')
setnames(tal.prior,'Plot','plotID')
setnames(tal.prior,'Perc.Soil.Moisture','moisture')
setnames(tal.prior,'Bulk.Density.','bulk.density')
tal.prior$soilTemp <- NA

#change horizon labels of AH/OH to M/O for talbot data to match neon.
tal.prior$horizon2 <- NA
tal.prior$horizon2 <- ifelse(tal.prior$horizon == 'AH','M',tal.prior$horizon2)
tal.prior$horizon2 <- ifelse(tal.prior$horizon == 'OH','O',tal.prior$horizon2)
tal.prior$horizon  <- tal.prior$horizon2
tal.prior <- tal.prior[,horizon2 := NULL]

#extract Ndep- note no coverage in Alaska.
tal.ndep <- extract_ndep(tal.prior$longitude,tal.prior$latitude)
tal.prior <- cbind(tal.prior, tal.ndep)

#get prism climate for talbot - note no coverage for alaska
#note - you need more prism files for older prism dates for this to work. 
#tal.climate <- prism_query(tal.prior, prism.dir)

#get PRISM climate data and ndep for NEON ITS
neon.climate <- prism_query(neon.map)
neon.ndep    <- extract_ndep(neon.map$longitude,neon.map$latitude)
neon.map     <- cbind(neon.map,neon.climate,neon.ndep)

#query soil physical properties to get horizon data.
soil_phys <- neon_query(neon.map,"DP1.10086.001")
neon.map <- merge(neon.map,soil_phys, all.x = T, by.x='geneticSampleID', by.y='geneticSampleID')

#get EM relative abundance for neon data.
neon.tax$ECM <- ifelse(neon.tax$guild == 'Ectomycorrhizal',1,0)
neon.tax$ECM <- ifelse(is.na(neon.tax$ECM),0,neon.tax$ECM)
neon.otu.ecm <- neon.otu
for(i in 1:ncol(neon.otu.ecm)){
  neon.otu.ecm[,i] <- neon.otu.ecm[,i]*neon.tax$ECM
}
#get sequence counts.
n.seqs <- colSums(neon.otu)
n.seqs.ecm <- colSums(neon.otu.ecm)
neon.map$relEM <- n.seqs.ecm / n.seqs

#get doy for neon data
neon.map$doy <- lubridate::yday(neon.map$date)

#get a site level data frame for site level variables - climate, etc.
tal.site <- tal.prior[,.(site,longitude,latitude,date,doy,mat30,map30,n.dep,dry.dep,wet.dep)]
tal.site <- tal.site[!duplicated(tal.site$site),]
neon.site <- data.table(neon.map)
neon.site <- neon.site[,.(site,longitude,latitude,date,doy,mat30,map30,n.dep,dry.dep,wet.dep)]
neon.site <- neon.site[!duplicated(neon.site$site),]

#subset these to be simpler for working tomorrow.
neon.map <- data.table(neon.map)
 tal.core <- as.data.frame(tal.prior[,.(site,plotID,horizon,relEM,relAM,relSAP,relWR,longitude,latitude,date,epoch_date,doy,mat30,map30,n.dep,dry.dep,wet.dep,pH,moisture,cn)])
neon.core <- as.data.frame(neon.map[,.(site,plotID,horizon,relEM,longitude,latitude,date,epoch_date,doy,mat30,map30,n.dep,dry.dep,wet.dep,soilTemp)])
 tal.site <- as.data.frame(tal.site)
neon.site <- as.data.frame(neon.site)

#save formatted outputs.
saveRDS( tal.core, tal.core.path)
saveRDS(neon.core,neon.core.path)
saveRDS( tal.site, tal.site.path)
saveRDS(neon.site,neon.site.path)