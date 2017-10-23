#This function takes a data frame that has site and year_month columns and a vector of NEON data products of interest.
#It returns a dataframe of unique site X year_month combinations from your data and whether or not there are data for it.


#' Checking neon product availability for a set of site by year_month combinations.
#'
#' @param some_data a data frame that includes the columns 'site' and 'year_month'
#' @param neon_products a vector of neon products to check data availability for.
#'
#' @return This function returns a data frame identifying which site by year_month combinations have data for a given NEON product.
#' @export
#'
#' @examples
#' map <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')
#' of_interest <- c('DP1.10078.001','DP1.10086.001')
#' check_neon_products(map,of_interest)

check_neon_products <- function(some_data,neon_products){
    #isolate site data on interest
    some_data <- data.table(some_data)
    sites <- some_data[,.(site,year_month)]
    sites$site_year_month <- paste(sites$site,sites$year_month,sep='_')
    setkey(sites,site_year_month)
    sites <- unique(sites)
    
    #NEON product loop
    product_list <- list()
    for(j in seq_along(neon_products)){
    this_neon_product <- neon_products[j]
    
        #site within product loop
        out.list <- list()
        for(i in seq_along(sites$site)){
          #specify which site and date we are talking about
          this_site       <- sites[i,site]
          this_year_month <- sites[i,year_month]
          
          #Connect to the NEON API, grab particular product.
          physical <- nneo::nneo_data(product_code = this_neon_product,
                                      site_code = this_site,
                                      year_month = this_year_month, package = "expanded")
          
          #Are there data? Its a T/F based on whether there are files listed or not.
          out <- nrow(physical$data$files)  > 0
          
          #save specific site by year-month data to output list.
          out.list[[i]] <- out
        }
        #merge list of T/F results into a single vector for each NEON product.
        product_list[[j]] <- do.call(rbind,out.list)
    }
    
    #merge together T/F vectors w/ site by year combinations, return output.
    product_list <- do.call(cbind,product_list)
    colnames(product_list) <- neon_products
    return(cbind(sites,product_list))
}