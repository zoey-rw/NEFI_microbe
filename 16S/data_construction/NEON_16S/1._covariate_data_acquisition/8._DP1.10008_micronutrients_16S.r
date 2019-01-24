# Get soil chemistry data and bulk density data for each NEON plot
# from initial characterization data - mostly 2016

rm(list=ls())
source('paths.r')
library(dplyr)

output.path <- dp1.10008_dp.10047_output_16S.path


#get site dates with 16S sequence data
site_dates <- readRDS(site_dates_16S.path)
site_output_chem <- list()
site_output_phys <- list()

#connect to NEON API for DP1.10008.00 - Soil chemical properties (Distributed initial characterization) 
req <- httr::GET(paste0("http://data.neonscience.org/api/v0/products/DP1.10008.001")) 
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates available for this product.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#subset prod.site_dates to the sites you actually care about for 16S sequence data.
prod.site_dates <- prod.site_dates[names(prod.site_dates) %in% names(site_dates)]

#loop through prod.site_dates, grabbing initial characterization chemistry data
for(i in 1:length(prod.site_dates)){
  site <- names(prod.site_dates)[i] #specify site.
  date_output <- list()
  
  for(k in 1:length(prod.site_dates[[i]])){  #loop through dates within a site.
    date <- prod.site_dates[[i]][k] #specify unique site-date being queried.
    site.date <- paste0(site,'/',date)
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
  
    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("spc_biogeochem", core.files$data$files$name)) == 0){
      next
    }
    core  <- read.delim(core.files$data$files$url
                        [intersect(grep("spc_biogeochem", core.files$data$files$name),
                                   grep("expanded", core.files$data$files$name))], sep=",")
    
    core$dateID <- date #write date to table
    date_output[[k]] <- core #save output in date list
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output_chem[[i]] <- date_output
}
site_output_chem <- do.call(plyr::rbind.fill, site_output_chem)


#connect to NEON API for DP1.10047.00 - Soil physical properties
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10047.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates availabile for this product.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#subset prod.site_dates to the sites you actually care about for 16S sequence data.
prod.site_dates <- prod.site_dates[names(prod.site_dates) %in% names(site_dates)]

#loop through prod.site_dates, grabbing initial characterization physical data
site_output <- list()
for(i in 1:length(prod.site_dates)){
  #specify site.
  site <- names(prod.site_dates)[i]
  date_output <- list()
  
  #loop through dates within a site.
  for(k in 1:length(prod.site_dates[[i]])){
    #specify unique site-date being queried.
    date <- prod.site_dates[[i]][k]
    site.date <- paste0(site,'/',date)
    
    #grab data for a particular site-date combination.
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))

    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("bulkdensity", core.files$data$files$name)) == 0){
      next
    }
    core  <- read.delim(core.files$data$files$url
                        [intersect(grep("bulkdensity", core.files$data$files$name),
                                   grep("basic", core.files$data$files$name))], sep=",")
    core$dateID <- date #write date to table
    date_output[[k]] <- core #save output in date list
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output_phys[[i]] <- date_output
}
site_output_phys <- do.call(plyr::rbind.fill, site_output_phys)

saveRDS(list(site_output_chem, site_output_phys), output.path)
