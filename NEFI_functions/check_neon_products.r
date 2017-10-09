#This function takes a data frame that has site and year_month columns and a vector of NEON data products of interest.
#It returns a dataframe of unique site X year_month combinations from your data and whether or not there are data for it.

#testing function.
#clear R environment, load packages.
#rm(list=ls())
#library(nneo)
#library(data.table)
#map <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')
#of_interest <- c('DP1.10078.001','DP1.10086.001')
#test <- check_neon_products(map,of_interest)
#head(test)

check_neon_products <- function(some_data,neon_products){
    #isolate site data on interest
    some_data <- data.table(some_data)
    sites <- some_data[,.(site,year_month)]
    sites$site_year_month <- paste(sites$site,sites$year_month,sep='_')
    setkey(sites,site_year_month)
    sites <- unique(sites)
    
    #NEON product loop
    product_list <- list()
    for(j in 1:length(neon_products)){
    this_neon_product <- neon_products[j]
    
        #site within product loop
        out.list <- list()
        for(i in 1:nrow(sites)){
          #specify white site and date we are talking about
          this_site       <- sites[i,site]
          this_year_month <- sites[i,year_month]
          
          #Connect to the NEON API, grab soil physical properties.
          physical <- nneo_data(product_code = this_neon_product,
                                site_code = this_site,
                                year_month = this_year_month, package = "expanded")
          
          #Are there data? Its a T/F based on whether there are files listed or not.
          if(nrow(physical$data$files)  > 0){out = T}
          if(nrow(physical$data$files) == 0){out = F}
          
          #save specific site by year-month data to output list.
          out.list[[i]] <- out
        }
        #merge list of T/F results into a single vector for each NEON product.
        product_list[[j]] <- do.call(rbind,out.list)
    }
    
    #merge together T/F vectors w/ site by year combinations, return output.
    product_list <- do.call(cbind,product_list)
    colnames(product_list) <- neon_products
    out <- cbind(sites,product_list)
    return(out)
}