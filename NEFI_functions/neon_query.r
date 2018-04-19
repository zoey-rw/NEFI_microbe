#' Query a neon product for a given set of site by year_month combinations, or all observations for a given set of sites.
#' Depends on nneo and data.table packages.
#' Function will handle redundant entries (i.e. multiple rows with the same site by year_month combination).
#' If you don't get data back for a particular site by year_month combination, then those data are not on the NEON server.
#' Currently depends on a grep_string lookup table. Hopefully this will be deprecated by a better solution.
<<<<<<< HEAD
#' Throws an error if you submit a data.frame that returns no data.
=======
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c
#'
#' @param some_data dataframe with the columns 'site' and 'year_month'
#' site should be formatted as a character entry in all caps, ex. "BART"
#' year_month should be formatted as YYYY-MM, ex. "2014-08"
#' 
#' @param neon_product_code product code for a given neon product
<<<<<<< HEAD
#' @param site_level logical, should this query be done at the site level? Default T. 
#' If you want a subset of dates, set to F and be sure to have a vector of year-months.
=======
#' @param site_level logical, should this query be done at the site level? Default F.
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c
#' If set to T, will return all observations made at any time point for a given data product at a given site.
#' @param lookup.path this links the product codes to grep strings for pulling out actual dataframe of interest from neon.
#'
#' @return returns a dataframe with all data present in a product for all site by year_month combinations.
#' @export
#'
#' @examples
#' map <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')
#' map <- map[c(1:5,(nrow(map)-5):nrow(map)),]
#' this.prod <- c("DP1.10086.001")
#' this.prod <- c("DP1.00097.001")
#' this.prod <- c("DP1.10078.001")
#' test <- neon_query(map, this.prod)
#' test <- neon_query(map, this.prod, site_level = T)

neon_query <- function(some_data,
                       neon_product_code,
<<<<<<< HEAD
                       site_level = T,
=======
                       site_level = F,
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c
                       lookup.path='/home/caverill/NEFI_microbe/neon_product_grep_lookup.csv'){
  
  #grab data frame that has only unique combinations of site and year_month.
  if(site_level == F){
    some_data <- data.table::data.table(some_data)
    sites <- some_data[,.(site,year_month)]
    sites$site_year_month <- paste(sites$site,sites$year_month,sep='_')
    data.table::setkey(sites,site_year_month)
    sites <- unique(sites)
  }
  
  #If you are pulling all data for a set of sites then we generate a different site by year_month table
  if(site_level == T){
    sites <- unique(some_data$site)
    all_dates <- list()
    for(i in seq_along(sites)){
    site_status <- nneo::nneo_site(sites[i])
    dates <- site_status$dataProducts[site_status$dataProducts$dataProductCode == neon_product_code,3]
    #This is a hack because some sites have an initial NA value in date list, others don't. 
    #I think this is a NEON bug.
    if(length(dates) > 1){
      dates[[1]] <- dates[[length(dates)]]
      dates <- dates[1]
    }
    if(length(dates) == 0){dates[[1]] <- numeric()}
    
    #if there are observations for certain dates put them together.
    if(length(dates[[1]]) > 0){
        dates_v <- do.call('rbind',dates)
        all_dates[[i]] <- data.frame(rep(sites[i],length(dates_v)),as.vector(dates_v))
        colnames(all_dates[[i]]) <- c('site','year_month')
    }
    #if there are not, don't.
    if(length(dates[[1]]) == 0){all_dates[[i]] = NULL}
    }
    
    #put it all together   
    sites <- do.call('rbind',all_dates)
    sites$site_year_month <- paste(sites$site,sites$year_month,sep='_')
    sites <- data.table::data.table(sites)
  }
  
  #grab the grep string for the particular product
  lookup <- read.csv(lookup.path, header=T)
  grep_string <- lookup[lookup$neon_product == neon_product_code,c('grep_string')]

  out.list <- list()
  for(i in seq_along(sites$site)){
    #specify white site and date we are talking about
    this_site       <- sites[i,site]
    this_year_month <- sites[i,year_month]
    
    #Connect to the NEON API, grab NEON product of interest for a site by year_month combination.
<<<<<<< HEAD
    site_product <- nneo::nneo_data(product_code = neon_product_code, 
                                            site = this_site, 
                                      year_month = this_year_month)
=======
    site_product <- nneo::nneo_data(product_code = neon_product_code, site = this_site, year_month = this_year_month)
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c

    #only download data if there are data. duh.
    if(nrow(site_product$data$files) > 0){
      #download the .csv file that actually has the data you want.
      path <- site_product$data$files$url[grep(grep_string,site_product$data$files$name)]
      #take last option in path list if more than one option. Some products have 'basic' and 'expanded' versions
      path <- path[length(path)]
<<<<<<< HEAD
      #currently running a try, as NEON sometimes throws "HTTP status was '500 Internal Server Error'" for some data products that used to work.
      out <- try(read.delim(path, sep = ','), silent = T)
    }
    
    #if there are no data, return NULL to list. If there is a server error, pass NULL as well.
    if(nrow(site_product$data$files) == 0){out = NULL}
    if(class(out) == 'try-error')         {out = NULL}
=======
      #take the last option in the path list. Some have basic and expanded versions.
      #currently running a try, as there are NEON is throwing "HTTP status was '500 Internal Server Error'" for some data products that used to work.
      #out <- read.delim(path, sep = ',')
      out <- try(read.delim(path, sep = ','), silent = T)
    }
    
    #if there are no data, return NA to list. If there is a server error, pass NULL as well.
    if(class(out) == 'try-error'){out = NULL}
    if(nrow(site_product$data$files) == 0){out = NULL}
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c
    
    #save specific site by year-month data to output list.
    out.list[[i]] <- out
  }
<<<<<<< HEAD

  #If there is output to return...
  if(length(out.list) > 0){
    #cull NA entries in output, merge together. for some reason a few obs 
    out.list <- out.list[!sapply(out.list,is.null)]
    #grab first object in out list. 
    z <- out.list[[1]]
    for(i in 2:length(out.list)){
      k <- out.list[[i]]
      #if  dataframe k doesn't have all the column names in dataframe z add them, set to NA.
      if(sum(colnames(k) %in% colnames(z)) < ncol(z)){
        to_add <- colnames(z[,!(colnames(z) %in% colnames(k))])
        k[,to_add] <- NA
      }
      #if data frame k has columns not in dataframe z, add then and set to NA.
      if(sum(colnames(z) %in% colnames(k)) < ncol(k)){
        to_add <- colnames(k[,!(colnames(k) %in% colnames(z))])
        z[,to_add] <- NA
      }
      #now that they def have same columns with same names, merge them. 
      z <- rbind(z,k)
    }    
  }
  
  #If there is not data to return, stop and say this.
  if(length(out.list) == 0){
    stop('No data to return for this product!')
  }

  #return output, which is a dataframe, or a message that says things broke.
=======
  #cull NULL entries in output, merge together. for some reason a few obs 
  out.list <- out.list[!sapply(out.list,is.null)]
  #out.list <- do.call('rbind',out.list)
  
  #try doing this another way, assigning NULL columns as necssary.
  z <- out.list[[1]]
  for(i in 2:length(out.list)){
    k <- out.list[[i]]
    #if the dataframe you are merging doesn't have all the column names so far add them, set to NULL.
    if(sum(colnames(k) %in% colnames(z)) < ncol(z)){
      to_add <- colnames(z[,!(colnames(z) %in% colnames(k))])
      k[,to_add] <- NA
    }
    #if data frame you are merging has columns not in zz, address this as well.
    if(sum(colnames(z) %in% colnames(k)) < ncol(k)){
      to_add <- colnames(k[,!(colnames(k) %in% colnames(z))])
      z[,to_add] <- NA
    }
    #now that they def have same columns with same names, merge them. 
    z <- rbind(z,k)
  }
  
  #return output, which is a dataframe.
>>>>>>> 86f6aad5e093b910fe7cc8ff293b9089204f114c
  return(z)
  
} #end function.