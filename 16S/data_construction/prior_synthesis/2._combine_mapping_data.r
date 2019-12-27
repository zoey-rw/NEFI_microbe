# combining Delgado and Ramirez covariate data
rm(list=ls())
source('NEFI_functions/crib_fun.r')
source("paths.r")
source("paths_fall2019.r")
source("NEFI_functions/convert_K.r")
source("NEFI_functions/extract_C.r")
library(runjags)

output.path <- delgado_ramirez_bahram_mapping.path

#### prep Ramirez et al. data ####
map.ram <- readRDS(ramirez_clean_map.path)
map.ram$source <- "Ramirez"
map.ram$study_id <- map.ram$dataset
#map.ram$pC <- extract_C(map.ram$latitude, map.ram$longitude, topsoil=T)/10
#map.ram$new.C.30 <- extract_C(map.ram$latitude, map.ram$longitude, topsoil=F)/10
map.ram <- map.ram[which(map.ram$sequencing_platform != "454"),]


# prep delgado mapping data
map.del <- readRDS(delgado_metadata_spatial.path)
map.del$source <- "Delgado"
map.del$study_id <- "Delgado"
map.del$depth_max <- 7.5
#map.del$new.C.5 <- extract_C(map.del$latitude, map.del$longitude, topsoil=T)/10
#map.del$new.C.30 <- extract_C(map.del$latitude, map.del$longitude, topsoil=F)/10


# fill with NAs? or, only common columns?
common_cols <- intersect(colnames(map.ram), colnames(map.del))
master.map <- do.call(plyr::rbind.fill, list(map.ram, map.del))
master.map <- master.map[,colnames(master.map) %in% c(common_cols)]

# add some colors
# study_id <- levels(as.factor(master.map$study_id))
# pals <- distinctColorPalette(length(study_id))
# colors <- cbind(study_id, pals)
# master.map <- merge(master.map, colors)

# subset to northern temperate latitudes
master.map <- master.map[master.map$latitude < 66.5 & master.map$latitude > 23.5,]

# reduce map values
master.map$map <- master.map$map/1000
master.map$map_sd <- master.map$map_sd/1000

# save output
saveRDS(master.map, output.path)
