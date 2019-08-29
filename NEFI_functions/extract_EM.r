# function for extracting EM % from maps from Steidinger et al. 2019
extract_EM <- function(latitude, longitude, path = "/projectnb/talbot-lab-data/NEFI_data/covariate_data/EM_AM_maps/em_1degreeMasked.nc") {
  # load EM proportion raster
  em_raster <- raster::raster(path)
  points <- cbind(longitude, latitude)
  # extract values using points
  output <- raster::extract(em_raster, points)
  return(output)
}