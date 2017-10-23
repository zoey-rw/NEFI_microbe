#get soil chemical data from NEON sites across all sampling dates.
#clear environment, load packages.
rm(list=ls())
library(nneo)

d <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')

test <- soil_chem_query(d)

soil_chem_query <- function(some_data){
    sites <- unique(d$site)
    prod_code <- 'DP1.00097.001'
    
    site_chem <- list()
    for(i in 1:length(sites)){
      #for each site grab all the dates that have soil chemical data.
      site_status <- nneo_site(sites[i])
            dates <- site_status$dataProducts[site_status$dataProducts$dataProductCode == prod_code,3]
            dates <- do.call('rbind',dates)
            
            #for each date grab the actual chemical data.
            site_date <- list()
            for(j in 1:length(dates)){
              chemical <- nneo_data(product_code = prod_code, site = sites[i], year_month = dates[j])
              if(nrow(chemical$data$files) > 0){
                  site_date[[j]]  <- read.delim(chemical$data$files$url
                                     [grep('mgc_perbiogeosample',chemical$data$files$name)], sep = ',')
              }
              #don't need to if there are no data for that date. Though this should never happen...
              if(nrow(chemical$data$files) == 0){site_date[[j]] = NULL}
            }
        #collapse output
        site_date <- do.call('rbind', site_date)
        site_chem[[i]] <- site_date
    }
    #rbind all the data across sites and rturn as output
    site_chem <- do.call('rbind', site_chem)
    return(site_chem)
} #end function