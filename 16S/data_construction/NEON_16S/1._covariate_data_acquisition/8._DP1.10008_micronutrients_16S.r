# Get micronutrient (Ca, Mg, P, K) data and bulk density data for each NEON plot

rm(list=ls())
source('paths.r')
library(dplyr)

#get site dates with 16S sequence data
site_dates <- readRDS(site_dates_16S.path)
site_output_chem <- list()
site_output_phys <- list()


#connect to NEON API for DP1.10008.00 - Soil chemical properties (Distributed initial characterization) 
req <- httr::GET(paste0("http://data.neonscience.org/api/v0/products/DP1.10008.001")) 
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates available for this product.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#subset prod.site_dates to the sites you actually care about for 16S sequence data.
prod.site_dates <- prod.site_dates[names(prod.site_dates) %in% names(site_dates)]

#loop through prod.site_dates, grabbing initial characterization chemistry data
for(i in 1:length(prod.site_dates)){
  site <- names(prod.site_dates)[i] #specify site.
  date_output <- list()
  
  for(k in 1:length(prod.site_dates[[i]])){  #loop through dates within a site.
    date <- prod.site_dates[[i]][k] #specify unique site-date being queried.
    site.date <- paste0(site,'/',date)
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))
  
    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("spc_biogeochem", core.files$data$files$name)) == 0){
      next
    }
    core  <- read.delim(core.files$data$files$url
                        [intersect(grep("spc_biogeochem", core.files$data$files$name),
                                   grep("expanded", core.files$data$files$name))], sep=",")
    
    core$dateID <- date #write date to table
    date_output[[k]] <- core #save output in date list
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output_chem[[i]] <- date_output
}
site_output_chem <- do.call(plyr::rbind.fill, site_output_chem)


#connect to NEON API for DP1.10047.00 - Soil physical properties
req <- httr::GET("http://data.neonscience.org/api/v0/products/DP1.10047.001")
req.text <- httr::content(req, as="text")
avail <- jsonlite::fromJSON(req.text, simplifyDataFrame=T, flatten=T)

#grab a vector of the urls to data. One per unique site-date combination.
urls <- unlist(avail$data$siteCodes$availableDataUrls)

#get a nested list of site by dates availabile for this product.
prod.site_dates <- avail$data$siteCodes[,2]
names(prod.site_dates) <- avail$data$siteCodes[,1]

#subset prod.site_dates to the sites you actually care about for 16S sequence data.
prod.site_dates <- prod.site_dates[names(prod.site_dates) %in% names(site_dates)]

#loop through prod.site_dates, grabbing initial characterization physical data
site_output <- list()
for(i in 1:length(prod.site_dates)){
  #specify site.
  site <- names(prod.site_dates)[i]
  date_output <- list()
  
  #loop through dates within a site.
  for(k in 1:length(prod.site_dates[[i]])){
    #specify unique site-date being queried.
    date <- prod.site_dates[[i]][k]
    site.date <- paste0(site,'/',date)
    
    #grab data for a particular site-date combination.
    core.JSON  <- httr::GET(urls[grep(site.date, urls)])
    core.files <- jsonlite::fromJSON(httr::content(core.JSON, as='text'))

    #check if there are even core data for a site-date combo. If not, skip!
    if(length(grep("bulkdensity", core.files$data$files$name)) == 0){
      next
    }
    core  <- read.delim(core.files$data$files$url
                        [intersect(grep("bulkdensity", core.files$data$files$name),
                                   grep("basic", core.files$data$files$name))], sep=",")
    core$dateID <- date #write date to table
    date_output[[k]] <- core #save output in date list
  }
  date_output <- do.call(plyr::rbind.fill, date_output)
  site_output_phys[[i]] <- date_output
}
site_output_phys <- do.call(plyr::rbind.fill, site_output_phys)

# merge nutrient data with bulk density data
merged <- merge(site_output_phys, site_output_chem, 
                by = c("horizonID"))

# drop unnecessary columns (see variables file for information on each)
nutr_df <- merged[,(colnames(merged) %in% c("uid.x", "siteID.x", "plotID.x", "collectDate.x", "horizonID",
                                                 "horizonName.x", "biogeoTopDepth", "biogeoBottomDepth",
                                                 "mgNh4d", "caNh4d", "kNh4d", "OlsenPExtractable", "pOxalate", 
                                                 "bulkDensTopDepth", "bulkDensBottomDepth", "bulkDensVolume", "bulkDensWetWeight",
                                                 "bulkDensDryWeight", "bulkDensOvenDry", "bulkDensFieldMoist"))]


# take the natural log of the nutrients
nutr_df[,c("LogCa","LogK","LogMg","LogP")] <- log(nutr_df[,c("caNh4d","kNh4d","mgNh4d","pOxalate")]) # or "OlsenPExtractable"

# calculate bulk density for the rows where it is empty
# dividing the dry weight (in grams) by the volume (in cm3)
# for some reason, rows that don't have calculated bulk density DO have volume + mass, and vice versa.
nutr_df$bd <- NA
nutr_df$bd <- ifelse(!is.na(nutr_df$bulkDensOvenDry), 
                     nutr_df$bulkDensOvenDry, 
                     nutr_df$bulkDensDryWeight/nutr_df$bulkDensVolume)

# take the horizons that do not go deeper than 30 inches.
nutr_df <- nutr_df[nutr_df$bulkDensBottomDepth <= 30,]
# or take the top two horizons:
# nutr_df <- tbl_df(nutr_df) %>% group_by(plotID.x) %>% top_n(-2, bulkDensTopDepth)

# get means of log-transformed nutrient values, weighted by bulk density.
nutr_tbl <- sapply(split(nutr_df, nutr_df$plotID.x), function(x) apply(x[, 21:24], 2, weighted.mean, 
                                          x$bd)) 
# check results dataframe
head(nutr_tbl) 

# alright, negative values, that's not good. our units must be different from Tedersoo's, but they don't state units in their methods.
# NEON nutrients are in centimolesPerKilogram, other than phosphorus, which is milligramsPerKilogram

# check against bahram data - these values are already log-transformed
# bahram <- readRDS(bahram_metadata.path)

# other method: calculating per-horizon stocks, summing them, dividing by total volume
# this won't work because we're missing so many volume measurements

# multiply concentration by bulk density
# nutr_df$Ca_by_bd <- nutr_df$LogCa * nutr_df$bd
# nutr_df$K_by_bd <- nutr_df$LogK * nutr_df$bd
# nutr_df$Mg_by_bd <- nutr_df$LogMg * nutr_df$bd
# nutr_df$P_by_bd <- nutr_df$LogP * nutr_df$bd

#sum the nutrient stocks from the top horizons, then divide by the sum of the volumes.
# nutr_tbl <- aggregate(cbind(nutr_df$Ca_by_bd, nutr_df$K_by_bd, nutr_df$Mg_by_bd, nutr_df$P_by_bd), 
#             by=list(Category=nutr_df$plotID.x), FUN=sum)
# colnames(nutr_tbl) <- c("plotID", "Ca", "K", "Mg", "P")
# nutr_tbl$volume <-  aggregate(nutr_df$bulkDensVolume, 
#                     by=list(Category=nutr_df$plotID.x), FUN=sum)

# compare the plot-level nutrient data we generated with the plots present in our microbial abundance data.
map <- readRDS(obs.table_16S.path) # read in observations
obs_plots <- unique(map$plotID) # get unique observation plotIDs
plots <- unique(nutr_tbl$plotID) # get unique nutrient plotIDs

# um. we are missing nutrient data for most of the plots that we have observations for.
setdiff(obs_plots, plots) # plotIDs only in our observations
setdiff(plots, obs_plots) # plotIDs only in the nutrient dataset
intersect(obs_plots, plots) # plotIDs common to both
# DAMN IT


# notes

# we are keeping data for: Ammonium acetate (AA) extractable magnesium, calcium, and potassium, and phosphorus
# but which phosphorus measurement we should actually keep? pOxalate, probably. look how many missing values there are for Olsen P.
# nrow(site_output_chem[is.na(site_output_chem$MehlichIIITotP),]) #128
# nrow(site_output_chem[is.na(site_output_chem$Bray1PExtractable),]) # 851
# nrow(site_output_chem[is.na(site_output_chem$OlsenPExtractable),]) # 443
# nrow(site_output_chem[is.na(site_output_chem$pOxalate),]) # 5

# If I'm interpreting correctly...
# P extracted with Ammonium oxalate is related to P from Ammonium lactate (priors) using this equation:
# Pox = 933 + 27 PAL
# From Table 5, https://skemman.is/bitstream/1946/19877/1/Gudmundsson%20et%20al%202014.pdf
# legitimate enough????

# AA-derived potassium measurements has a .97 correlation with Ammonium lactate (AL)-derived measurements from priors,
# Zebec et al. 2017: https://link.springer.com/content/pdf/10.1134%2FS1064229317130051.pdf

# Olsen phosphorus has a .908 correlation with our prior measurements of AL-derived phosphorus: 
# Horta et al 2010, https://www.tandfonline.com/doi/abs/10.1080/00103624.2010.508296 

# which bulk density columns should we keep?
# bulkDensVolume,Volume of the bulk density sample,real,cubicCentimeter,basic
# bulkDensWetWeight,Total weight of the bulk density sample prior to drying,real,gram,basic
# bulkDensDryWeight,Total weight of the bulk density sample after drying,real,gram,basic
# bulkDensCoarseFragWeight,Weight of the coarse (>2 mm) fragments in the bulk density sample,real,gram,basic
# bulkDensCoarseFragDens,Density of the coarse (>2 mm) fragments in the bulk density sample,real,gramsPerCubicCentimeter,basic
# bulkDensThirdBar,Weight per unit volume measured after equilibration at one third bar water tension,real,gramsPerCubicCentimeter,basic
# bulkDensOvenDry,Weight per unit volume measured on oven dried soil clods,real,gramsPerCubicCentimeter,basic
# bulkDensFieldMoist,Weight per unit volume measured at field moisture,real,gramsPerCubicCentimeter,basic
# waterRetentionThirdBar,Water content after equilibration at one-third bar water tension, reported as gravimetric percent on the <2 mm fraction,real,percent,basic
# fieldWaterContent,Water content under field conditions, reported as gravimetric percent on the <2 mm fraction,real,percent,basic
