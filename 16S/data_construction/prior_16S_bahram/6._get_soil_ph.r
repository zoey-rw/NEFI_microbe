# get KCl and CaCl pH values from the World Soil Information Service (WoSIS) 
# for converting NEON pH values to match prior values

rm(list=ls())
library(rgdal)
library(gdalUtils)

# downloading these using a "web feature service" (WFS) connection
dsn <- "WFS:http://data.isric.org/geoserver/wosis_latest/wfs/"

# read in KCl pH data
ogr2ogr(dsn, "wosis_latest_phkc", "wosis_latest:wosis_latest_phkc")
phkc <- readOGR("wosis_latest_phkc", "wosis_latest:wosis_latest_phkc")
kcl <- as.data.table(phkc@data)
kcl_alone <- kcl[,c("profile_id","phkc_val_1","layer_name", "lower_dept")]

# read in CaCl pH data
ogr2ogr(dsn, "wosis_latest_phca", "wosis_latest:wosis_latest_phca")
ph_cacl <- readOGR("wosis_latest_phca", "wosis_latest:wosis_latest_phca")
cacl <- as.data.table(ph_cacl@data)
cacl_alone <- cacl[,c("profile_id","phca_val_1","layer_name", "lower_dept")]

merged <- merge(cacl_alone, kcl_alone, by=c("profile_id", "layer_name", "lower_dept"))
merged <- merged[which(merged$lower_dept < 30),]
# dim(merged) # n = 7485

# create linear model of these two pH values
linearMod <- lm(phca_val_1 ~ phkc_val_1, data=merged)
print(linearMod) 
# consistent with other values in the literature, such as
# https://esdac.jrc.ec.europa.eu/Library/Data/PH/Documents/pH_Pub.pdf 
# and in DOIs: 10.1016/j.geodrs.2018.e00185 and 10.4141/cjss81-067 

# just checking that it looks ok
plot(merged$phkc_val_1, merged$phca_val_1, col="red")
abline(linearMod) 

# save dataframe of values
saveRDS(merged, "pH_conversion_data.rds")