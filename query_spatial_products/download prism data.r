#downloading PRISM climate data.
#use prism package for R.
rm(list=ls())
library(prism)

#set prism directory. This is currently in colin's data directory, will change.
prism.path <- "/fs/data3/caverill/PRISM"
options(prism.path = prism.path)

#get normals
get_prism_normals(type="tmean",resolution = "800m",annual=T, keepZip=F)
get_prism_normals(type=  "ppt",resolution = "800m",annual=T, keepZip=F)

#get all 2013-2016 data: tmean, tmin, tmax, ppt.
get_prism_monthlys(type = 'tmean', year = c(2013:2016), mon = c(1:12), keepZip=F)
get_prism_monthlys(type =   'ppt', year = c(2013:2016), mon = c(1:12), keepZip=F)
get_prism_monthlys(type =  'tmin', year = c(2013:2016), mon = c(1:12), keepZip=F)
get_prism_monthlys(type =  'tmax', year = c(2013:2016), mon = c(1:12), keepZip=F)