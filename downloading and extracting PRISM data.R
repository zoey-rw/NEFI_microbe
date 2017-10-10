#writing function to get mean annual temperature and precipitation from PRISM.
#requires some packages
rm(list=ls())
library(sp)
library(raster)
library(data.table)

#grab PRISM rasters
map.PRISM <- raster('/fs/data3/caverill/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil')
mat.PRISM <- raster('/fs/data3/caverill/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil')

#grab your data.
d <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')



#try with prism R package
library(prism)
options(prism.path = "/fs/data3/caverill/PRISM")

#get normals
get_prism_normals(type="tmean",resolution = "800m",annual=T, keepZip=F)
get_prism_normals(type="ppt" ,resolution = "800m",annual=T, keepZip=F)

#get unique year-month combinations
dates <- data.frame(unique(d$year_month))
colnames(dates) <- c('year_month')
dates <- data.table(dates)
dates$year  <- substr(dates$year_month,1,4)
dates$month <- substr(dates$year_month,6,7)

#download year-month data you actually need
for(i in 1:length(unique(dates$year))){
  this.year <- as.numeric(unique(dates$year)[i])
  months <- as.numeric(unique(dates[year == this.year,month]))
  for(j in 1:length(unique(months))){
      get_prism_monthlys(type='tmean', year=this.year, mon = months[j], keepZip=F)
     #get_prism_monthlys(type='ppt'  , year=this.year, mon = months[j], keepZip=F)
  }
}




test <- prism_query(d)
prism_query <- function(data){
  #place paths of raster files.
  map.PRISM <- raster('/fs/data3/caverill/PRISM_ppt_30yr_normal_800mM2_annual_bil.bil')
  mat.PRISM <- raster('/fs/data3/caverill/PRISM_tmean_30yr_normal_800mM2_annual_bil.bil')
  
  #grab longitude and latitude from data frame.
  points <- cbind(data$longitude, data$latitude)
  
  #extract products, append to data frame.
  map_out <- extract(map.PRISM,points)
  mat_out <- extract(mat.PRISM,points)
  
  #return a dataframe with the climate variables
  output <- data.frame(map_out,mat_out)
  colnames(output) <- c('MAP','MAT')
  return(output)
  }
