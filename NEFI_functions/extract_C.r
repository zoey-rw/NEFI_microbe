#' extract_C.r
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
#' @return           #returns a vector of C values.
#' @export
#'
#' @examples
#' 4 points from very different places as (longitude, latitude).
#' points <- structure(c(-102.644235, -58.015319, 50.508187, -148.2747566, 38.816171, -7.719157, 21.754477, 64.766184), .Dim = c(4L,2L))
#' test.out <- extract_C(points[,2], points[,1])
#' 
extract_C <- function(latitude, longitude,folder = '/projectnb/talbot-lab-data/NEFI_data/covariate_data/soil_c/',topsoil=FALSE){
  #default directory path is scc1. Check if you are pecan2. Can be modified to work on more local machines this way.
  host <- system('hostname', intern=T)
  if(host == 'pecan2'){folder <- '/fs/data3/caverill/'}
  
  #specify raster path
  
  raster.path <- ifelse(topsoil==TRUE, 
                        paste0(folder,'ORCDRC_M_sl2_250m_ll.tif'), 
                        paste0(folder,'ORCDRC_M_sl4_250m_ll.tif'))
  
  #get lat-long as an object.
  points <- cbind(longitude, latitude)
  
  #load SoilGrids 250m soil C raster.
  dat <- raster::raster(raster.path)
  
  #extract soil C values.
  output <- raster::extract(dat, points)
  
  #return output
  return(output)
}
