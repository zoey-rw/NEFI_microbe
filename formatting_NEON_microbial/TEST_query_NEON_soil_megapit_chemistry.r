#get site-level soil chemical data from NEON sites across all sampling dates.
#clear environment, load packages.
rm(list=ls())
source('NEFI_functions/neon_query.r')

d <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')

test <- neon_query(d,neon_product_code = c("DP1.10086.001"), site_level = T)
test <- neon_query(d,neon_product_code = "DP1.00097.001", site_level = F)

test <- soil_chem_query(d, "DP1.00097.001")


#old function that successfully queried soil megapit sampling. 
soil_chem_query <- function(some_data,
                            neon_product_code,
                            lookup.path='/home/caverill/NEFI_microbe/neon_product_grep_lookup.csv'){
    #grab unique sites
    sites <- unique(some_data$site)
    
    #grab the grep string for the particular product
    lookup <- read.csv(lookup.path, header=T)
    grep_string <- lookup[lookup$neon_product == neon_product_code,c('grep_string')]
    
    site_data <- list()
    for(i in 1:length(sites)){
      #for each site grab all the dates that have data for a given product.
      site_status <- nneo::nneo_site(sites[i])
            dates <- site_status$dataProducts[site_status$dataProducts$dataProductCode == neon_product_code,3]
            dates <- do.call('rbind',dates)
            
            #for each date grab the actual chemical data.
            site_date <- list()
            for(j in 1:length(dates)){
              chemical <- nneo::nneo_data(product_code = neon_product_code, site = sites[i], year_month = dates[j])
              if(nrow(chemical$data$files) > 0){
                  site_date[[j]]  <- read.delim(chemical$data$files$url
                                     [grep(grep_string,chemical$data$files$name)], sep = ',')
              }
              #don't need to if there are no data for that date. Though this should never happen...
              if(nrow(chemical$data$files) == 0){site_date[[j]] = NULL}
            }
        #collapse output
        site_date <- do.call('rbind', site_date)
        site_data[[i]] <- site_date
    }
    #rbind all the data across sites and turn as output
    site_data <- do.call('rbind', site_data)
    return(site_data)
} #end function