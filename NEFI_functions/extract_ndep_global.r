#' Extract wet + dry N deposition data from Ackerman et al. 
#' Reported in kg per hectare per year, for the year 2014.
#' Not vectorized
#'
#' @param longitude a vector of site longitude
#' @param latitude  a vector of site latitude
#' @param folder path to a CSV file of annual inorganic nitrogen per km
#' folder path currently defaults to the directory in the lab's folder on the SCC
#'
#' @return returns a list of total ndeposition values for 2014 for all locations.
#' @export
#' @import 
#'
#'# from Ackerman et al.: https://conservancy.umn.edu/handle/11299/197613
extract.ndep.global <- function(latitude, longitude, 
                                folder = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/inorganic_N_deposition.csv"){
  require(MALDIquant)
ndep.glob <- read.csv(folder)
ndep.tot.2014 <- ndep.glob[,c("latitude", "longitude", "pixel_area_km2","tot_2014")]
ndep.global <- list()
for (n in 1:length(latitude)){
  lat <- latitude[n]
  lon <- longitude[n]
  sorted_lat <- sort(ndep.tot.2014$latitude)
  newlat <- sorted_lat[MALDIquant::match.closest(lat, sorted_lat)]
  sorted_lon <- sort(ndep.tot.2014$longitude)
  newlon <- sorted_lon[MALDIquant::match.closest(lon, sorted_lon)]
  ndep.point <- ndep.tot.2014[which(ndep.tot.2014$latitude==newlat & ndep.tot.2014$longitude==newlon),]
  ndep.global[n] <- ndep.point$tot_2014/100 #convert per-km to per-hectare
}  
output <- unlist(ndep.global)
return(output)
}
