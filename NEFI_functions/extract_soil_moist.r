#' extract_soil_moist.r
#' Extracts soil moisture from the Atlas of the Biosphere product, a MODIS NPP product.
#' Willmott, Cort J., and Matsuura Kenji (2001): Terrestrial Water Budget Data Archive: Monthly Time Series (1950-1999). 
#'
#' @param latitude   #vector of latitude in decimal degrees.
#' @param longitude  #vector of longitude in decimal degrees.
#' @param path       #path to soil moisture raster with header file. downloaded here: https://nelson.wisc.edu/sage/data-and-models/atlas/data/soilmoisture.zip
#' @return           #returns a vector of soil moisture values.
#' @export
#'
#' @examples
#' points <- structure(c(-102.644235, -58.015319, 50.508187, -148.2747566, 38.816171, -7.719157, 21.754477, 64.766184), .Dim = c(4L,2L))
#' test.out <- extract_soil_moist(points[,2], points[,1])
extract_soil_moist <- function(latitude, longitude, path = '/fs/data3/caverill/AoB_NPP/soil_moisture_raster/soilmoisture/w001001.adf'){
  #load soil moisture raster
  moist_raster <- raster::raster(path)
  #extract
  points <- cbind(longitude, latitude)
  output <- raster::extract(moist_raster, points)
  return(output)
}