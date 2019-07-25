#subset mapping file and otu table to only include soil samples from 'natural' ecosystems.
#clear environmenta, load packages.
rm(list=ls())
source('paths.r')
library(data.table)

raw.map <- fread(emp_map.path)

#Subset to soil observations.
map <- raw.map[empo_3 == 'Soil (non-saline)' & envo_biome_1 == 'terrestrial biome',]

#we probably want to make some more exclusion choices. There are agricultural sites here for instance.
#Lets exclude crops, polar desert, rangeland and "".
nope <- c('cropland biome','','polar desert biome','rangeland biome')
map <- map[!(envo_biome_3 %in% nope),]

#exclude a few studies that are sand from stream water filters and agroforestry.
nope <- c(755,1711,1714)
map <- map[!(study_id %in% nope),]
map$sampleID <- map[,"#SampleID"]
#map <- map[!grep("slow_sand filter|agriculture|Agricultural", paste(map$Description,map$title)),]
setnames(map, c("latitude_deg", "longitude_deg"), c("latitude", "longitude"))

map <- as.data.frame(map)
#save mapping file.
saveRDS(map, emp_map_clean.path)


soils <- map
mapWorld <-
  borders("world", colour = "gray50", fill = "gray50") # create a layer of borders
soils.map <- ggplot() +
  ggtitle("4553 samples from Soil and terrestrial biome") +
  mapWorld +
  geom_point(aes(
    x = as.numeric(soils$longitude_deg),
    y = as.numeric(soils$latitude_deg)
  ),
  color = "blue",
  size = 3)
soils.map


nope <- c('cropland biome','','polar desert biome','rangeland biome')
non.ag <- map[!(envo_biome_3 %in% nope),]
non.ag.map <- ggplot() +
  ggtitle("2356 samples from Soil and terrestrial biome \n+ excluded cropland/polar desert/rangeland biomes)") +
  mapWorld +
  geom_point(aes(
    x = as.numeric(non.ag$longitude_deg),
    y = as.numeric(non.ag$latitude_deg)
  ),
  color = "blue",
  size = 3)
non.ag.map


#exclude a few studies that are sand from stream water filters and agroforestry.
nope <- c(755,1711,1714)
non.sf <- non.ag[!(study_id %in% nope),]
non.sf.map <- ggplot() +
  ggtitle("1386 samples from Soil and terrestrial biome \n+ excluded cropland/polar desert/rangeland biomes \n+ excluded studies with samples from stream water filters or agroforestry") +
  mapWorld +
  geom_point(aes(
    x = as.numeric(non.sf$longitude_deg),
    y = as.numeric(non.sf$latitude_deg)
  ),
  color = "blue",
  size = 3)
non.sf.map



esv.sub <- readRDS(emp_map_clean.path)
esv.map <- ggplot() +
  ggtitle("1130 samples that overlap with sequence data") +
  mapWorld +
  geom_point(aes(
    x = as.numeric(esv.sub$longitude_deg),
    y = as.numeric(esv.sub$latitude_deg)
  ),
  color = "blue",
  size = 3)
esv.map

library(gridExtra)
grid.arrange(soils.map, non.ag.map, non.sf.map, esv.map)


esv.ph <- esv.sub[!is.na(esv.sub$ph),]
esv.map <- ggplot() +
  ggtitle("1130 samples that overlap with sequence data") +
  mapWorld +
  geom_point(aes(
    x = as.numeric(esv.sub$longitude_deg),
    y = as.numeric(esv.sub$latitude_deg)
  ),
  color = "blue",
  size = 3)
esv.map
