#This function takes a data frame that has the columns "site" and "year_month". Can have other columns, they will be ignored.
#site should be formatted as a character entry in all caps, ex. "BART"
#year_month should be formatted as YYYY-MM, ex. "2014-08"
#Function will handle redundant entries (i.e. multiple rows with the same site by year_month combination).
#Will return a data.table of all variables Colin Averill thinks are interesting from this data product.
#If you don't get data back for a particular site by year_month combination, then those data are not on the NEON server.

soil_physical_query <- function(some_data){
  #grab data frame that has only unique combinations of site and year_month.
  some_data <- data.table(some_data)
  sites <- some_data[,.(site,year_month)]
  sites$site_year_month <- paste(sites$site,sites$year_month,sep='_')
  setkey(sites,site_year_month)
  sites <- unique(sites)
  
  out.list <- list()
      for(i in 1:nrow(sites)){
        #specify white site and date we are talking about
        this_site       <- sites[i,site]
        this_year_month <- sites[i,year_month]
        
        #Connect to the NEON API, grab soil physical properties.
        physical <- nneo_data(product_code = "DP1.10086.001",
                              site_code = this_site,
                              year_month = this_year_month, package = "expanded")
        
        #only download data if there are data. duh.
        if(nrow(physical$data$files) > 0){
          #download the .csv file that actually has the data you want.
          prod <- read.delim(physical$data$files$url
                             [grep('soilCore',physical$data$files$name)], sep = ',')
          prod <- data.table(prod)
          
          #subset information you actually want from this table.
          out <- prod[,.(sampleID,geneticSampleID,siteID,plotID, #IDs
                         horizon,soilTemp,litterDepth,sampleTopDepth,sampleBottomDepth, #soil information
                         coreCoordinateX,coreCoordinateY,decimalLatitude,decimalLongitude,elevation, #position information
                         collectDate)] #temporal information
        }
        
        #if there are no data, return NA to list
        if(nrow(physical$data$files) == 0){out = NULL}
        
        #save specific site by year-month data to output list.
        out.list[[i]] <- out
      }
  #cull NULL entries in output, merge together
  out.list <- out.list[!sapply(out.list,is.null)]
  out.list <- do.call('rbind',out.list)
  
  #return output list, which is now a dataframe, and a data.table
  return(out.list)
} #end function.