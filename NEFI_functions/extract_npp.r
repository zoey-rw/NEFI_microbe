#' extract_npp.r
#' Extracts NPP from the Atlas of the Biosphere product, a MODIS NPP product.
#' Foley, J.A., I.C. Prentice, N. Ramankutty, S. Levis, D. Pollard, S. Sitch, and A. Haxeltine (1996) An Integrated Biosphere Model of Land Surface Processes, Terrestrial Carbon Balance and Vegetation Dynamics, Global Biogeochemical Cycles, 10, 603-628.
#'
#' @param latitude   #vector of latitude in decimal degrees.
#' @param longitude  #vector of longitude in decimal degrees.
#' @param path       #path to NPP raster with header file. downloaded here: https://nelson.wisc.edu/sage/data-and-models/atlas/data/npp.zip
#' @return           #returns a vector of NPP values.
#' @export
#'
#' @examples
#' points <- structure(c(-102.644235, -58.015319, 50.508187, -148.2747566, 38.816171, -7.719157, 21.754477, 64.766184), .Dim = c(4L,2L))
#' test.out <- extract_npp(points[,2], points[,1])
extract_npp <- function(latitude, longitude, path = '/fs/data3/caverill/AoB_NPP/npp/w001001.adf'){
  #load NPP raster
  npp_raster <- raster::raster(path)
  #extract
  points <- cbind(longitude, latitude)
  output <- raster::extract(npp_raster, points)
  return(output)
}