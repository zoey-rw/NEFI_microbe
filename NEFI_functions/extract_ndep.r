#This function takes longitude and latitude (in that order) as inputs.
#returns a dataframe of 2000-2015 averages on wet N deposition, dry deposition, and total N deposition.
#There are currently warnings. ignore them. 

extract_ndep <- function(longitude,latitude){
  #load dry deposition rasters
  dry.list <- list()
  for(i in 0:15){
    dry.list[[i+1]] <- raster(paste0('/fs/data3/caverill/CASTNET_Ndep/dry_dep/n_dw-',2000 + i,'.e00'))
  }
  #load wet deposition rasters
  wet.list <- list()
  for(i in 0:15){
    wet.list[[i+1]] <- raster(paste0('/fs/data3/caverill/CASTNET_Ndep/wet_dep/n_ww-',2000 + i,'.e00'))
  }
  
  #Fix the crs of each raster.
  for(i in 1:length(dry.list)){
    crs(dry.list[[i]]) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23.0 +lon_0=-96.0 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }
  for(i in 1:length(wet.list)){
    crs(wet.list[[i]]) <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23.0 +lon_0=-96.0 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
  }
  
  #grab long/lat. (important longitude before latitude)
  points <- cbind(longitude,latitude)
  #reproject points
  points <- SpatialPoints(points, proj4string = CRS("+init=epsg:4326"))
  
  #extract dry deposition, convert to mean annual
  dry.out <- list()
  for(i in 1:length(dry.list)){
    dry.out[[i]] <- extract(dry.list[[i]], points)
  }
  dry.dep <- Reduce('+',dry.out) / length(dry.out)
  #extract wet deposition
  wet.out <- list()
  for(i in 1:length(wet.list)){
    wet.out[[i]] <- extract(wet.list[[i]], points)
  }
  wet.dep <- Reduce('+',wet.out) / length(wet.out)
  
  #total ndep and output.
  n.dep <- dry.dep + wet.dep
  output <- data.frame(n.dep, dry.dep, wet.dep)
  
  return(output)
}