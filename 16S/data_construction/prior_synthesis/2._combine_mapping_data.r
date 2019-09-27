# combining EMP, Ramirez, and Bahram covariate data
rm(list=ls())
source('NEFI_functions/crib_fun.r')
source("paths.r")
source("paths_fall2019.r")
source("NEFI_functions/convert_K.r")
source("NEFI_functions/extract_C.r")
library(randomcoloR)
library(runjags)

output.path <- delgado_ramirez_bahram_mapping.path

#### prep Ramirez et al. data ####
map.ram <- readRDS(ramirez_clean_map.path)
map.ram$source <- "Ramirez"
map.ram$study_id <- map.ram$dataset
map.ram$new.C.5 <- extract_C(map.ram$latitude, map.ram$longitude, topsoil=T)/10
map.ram$new.C.30 <- extract_C(map.ram$latitude, map.ram$longitude, topsoil=F)/10

#### prep Bahram data ####
map.bahram <- data.table(readRDS("/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/old/prior_abundance_mapping/Bahram/bahram_metadata.rds"))
map.bahram$source <- "Bahram"
map.bahram$sampleID <- map.bahram$Run
map.bahram$new.C.5 <- extract_C(map.bahram$latitude, map.bahram$longitude, topsoil=T)/10
map.bahram$new.C.30 <- extract_C(map.bahram$latitude, map.bahram$longitude, topsoil=F)/10

# prep delgado mapping data
map.del <- readRDS(delgado_metadata_spatial.path)
map.del$source <- "Delgado"
map.del$study_id <- "Delgado"
map.del$new.C.5 <- extract_C(map.del$latitude, map.del$longitude, topsoil=T)/10
map.del$new.C.30 <- extract_C(map.del$latitude, map.del$longitude, topsoil=F)/10
map.del$ph <- map.del$pH
#map.bahram$ph <- map.bahram$pH

map.bahram$study_id <- factor("Bah")
pH_conversion <- readRDS(pH_conversion.path)
ph_dat <- as.matrix(pH_conversion[,c("phkc_val_1","phaq_val_1")])
converted.ph <- convert_K(map.bahram$ph, ph_dat)
map.bahram$ph <- converted.ph$mean
# 
# biome_key <- data.frame(emp = c("forest biome", "grassland biome", "shrubland biome", NA, NA, NA, NA, NA, "anthropogenic terrestrial biome", NA, NA, NA),
#                         ram = c("forest", "grassland", "xeric shrubland", "peatland", "glacier forefield", "tundra", NA, NA, NA, NA, NA, NA),
#                         bah = c("Boreal forests","Grassland and shrubland", NA, NA, NA, NA, NA, "Dry tropical forests", "Mediterranean", "Moist tropical forests", "Temperate deciduous forests", "Temperate coniferous forests"),
#                         biome = c("forest", "grassland_shrubland", "grassland_shrubland", "other","other","other","other","forest","other","forest","forest","forest"))
# 
# map.bahram$biome <- biome_key$biome[match(map.bahram$Biome, biome_key$bah)]
# map.emp$biome <- biome_key$biome[match(map.emp$envo_biome_2, biome_key$emp)]
# map.ram$biome <- biome_key$biome[match(map.ram$habitat, biome_key$ram)]
#map.bahram <- map.bahram[,.(sampleID,ph,NPP,map,mat,aridity,longitude,latitude,mdr,source)]

# fix some column names for the sake of merging
# setnames(map.ram, old = c("altitude.elevation..meter.", "dataset"), new = c("elevation_m","study_id"))
# setnames(map.bahram, old = c("P", "pC", "pN"), new = c("p","c","n"))
# map.ram$doy <- strftime(map.ram$collection_date, format = "%j")
# map.emp$collection_timestamp <- substr(map.emp$collection_timestamp, 1, 10)
# map.emp$collection_timestamp[map.emp$collection_timestamp==""] <- NA
# map.emp$doy <- strftime(map.emp$collection_timestamp, format = "%j")
#map.emp$c <- NA
# to.merge_map.emp <- map.emp[,colnames(map.emp) %in% c("sampleID","study_id","source","n.dep","map","mat","aridity","mdr","ph","ph_estim","elevation_m","env_biome","env_feature","envo_biome_2","latitude_deg","longitude_deg","NPP")]
# to.merge_map.emp$ph_estim.1 <- NULL
# to.merge_map.ram <- map.ram[,colnames(map.ram) %in% c("sampleID", "dataset", "source","n.dep","map","mat","aridity","mdr","ph","ph_estim","altitude.elevation..meter.","n","c",'current_vegetation', 'habitat', 'latitude','longitude',"NPP")]

# fill with NAs? or, only common columns?
common_cols <- intersect(colnames(map.ram), colnames(map.bahram))
common_cols <- intersect(common_cols, colnames(map.del))
master.map <- do.call(plyr::rbind.fill, list(map.ram, map.del, map.bahram))
master.map <- master.map[,colnames(master.map) %in% c(common_cols)]

# create new pH object, combining measured with estimated
master.map$new.ph <- master.map$ph
master.map[is.na(master.map$new.ph),]$new.ph <- master.map[is.na(master.map$new.ph),]$ph_estim

# add some colors
study_id <- levels(as.factor(master.map$study_id))
pals <- distinctColorPalette(length(study_id))
colors <- cbind(study_id, pals)
master.map <- merge(master.map, colors)

# save output
saveRDS(master.map, output.path)
