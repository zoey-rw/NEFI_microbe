#' extract_pH.r
#' Supplied with latitude and longtiude, this returns pH values at 250m resolution from the SoilGrids data product.
#' http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0169748
#' https://soilgrids.org/
#' 
#' Depends on the following packages:
#' raster,
#'
#' @param latitude   #a vector of latitude
#' @param longitude  #a vector of longitude
#' @param folder     #path to folder with SoilGrids raster. Defaults to scc1, built in check to work on pecan2.
#'
#' @return           #returns a vector of pH values.
#' @export
#'
#' @examples
#' 4 points from very different places as (longitude, latitude).
#' points <- structure(c(-102.644235, -58.015319, 50.508187, -148.2747566, 38.816171, -7.719157, 21.754477, 64.766184), .Dim = c(4L,2L))
#' test.out <- extract_pH(points[,2], points[,1])
extract_pH <- function(latitude, longitude,folder = '/project/talbot-lab-data/spatial_raster_data/SoilGrids_uncertainty'){
  #default directory path is scc1. Check if you are pecan2. Can be modified to work on more local machines this way.
  host <- system('hostname', intern=T)
  if(host == 'pecan2'){folder <- '/fs/data3/caverill/NEFI_microbial_data/SoilGrids_uncertainty/'}
  
  #specify raster path
  raster.path <- paste0(folder,'PHIHOX_M_sl2_250m.tif')
  
  #get lat-long as an object.
  points <- cbind(longitude, latitude)
  
  #load SoilGrids 1km pH raster.
  dat <- raster::raster(raster.path)
  
  #extract pH values.
  output <- raster::extract(dat, points)
  
  #divide by ten because to transform to true pH.
  output <- output/10
  
  #return output
  return(output)
}
