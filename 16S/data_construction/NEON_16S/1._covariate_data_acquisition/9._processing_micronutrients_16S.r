# Get micronutrient (Ca, Mg, P, K) data for each NEON plot
# from initial characterization data - mostly 2016
# NEON data is processed to better match the methods of the priors
# Ca and Mg - same methods, only units change to match priors.
# K - method changes and units change.
# P - methods change, units stay the same.
# P and K will each have a mean, uncertainty, and variance.

rm(list=ls())
source('paths.r')
source("NEFI_functions/convert_K.r")
source("NEFI_functions/convert_P.r")
library(dplyr)

# Read in soil chemistry and soil physical data
site_output_chem <- readRDS(dp1.10008_dp.10047_output_16S.path)[[1]]
site_output_phys <- readRDS(dp1.10008_dp.10047_output_16S.path)[[2]]

# Read in raw data comparing P-extraction methods
# oxalate-olsen data from Wuenscher: https://www.agriculturejournals.cz/publicFiles/144700.pdf
oxalate_olsen <- read.csv(paste0(pecan_gen_16S_dir, "P_oxalate_olsen_Wunscher.csv"))
# olsen-AL data from Erikkson: https://stud.epsilon.slu.se/826/1/eriksson_a_k_100209.pdf
olsen_AL <- read.csv(paste0(pecan_gen_16S_dir, "P_AL_olsen_mehlich_Erikkson.csv"))


#### Merge nutrient data with bulk density data ####

merged <- merge(site_output_phys, site_output_chem, 
                by = c("horizonID"))

# Drop unnecessary columns (see NEON variables files for information on each)
nutr_df <- merged[,(colnames(merged) %in% c("uid.x", "siteID.x", "plotID.x", "collectDate.x", "horizonID",
                                            "horizonName.x", "biogeoTopDepth", "biogeoBottomDepth",
                                            "mgNh4d", "caNh4d", "kNh4d", "OlsenPExtractable", "pOxalate", 
                                            "bulkDensTopDepth", "bulkDensBottomDepth", "bulkDensVolume", "bulkDensWetWeight",
                                            "bulkDensDryWeight", "bulkDensOvenDry", "bulkDensFieldMoist"))]

# Calculate bulk density for the rows where it is empty
# dividing the dry weight (in grams) by the volume (in cm3)
# rows that don't have calculated bulk density DO have volume + mass, and vice versa. Zoey doesn't know why.
nutr_df$bd <- NA
nutr_df$bd <- ifelse(!is.na(nutr_df$bulkDensOvenDry), 
                     nutr_df$bulkDensOvenDry, 
                     nutr_df$bulkDensDryWeight/nutr_df$bulkDensVolume)

# subset to the horizons that do not go deeper than 30 inches.
nutr_df <- nutr_df[nutr_df$bulkDensBottomDepth <= 30,]
# or take the top two horizons:
# nutr_df <- tbl_df(nutr_df) %>% group_by(plotID.x) %>% top_n(-2, bulkDensTopDepth)

#### Convert phosphorus ####

# subset conversion data to necessary columns
oxalate_olsen <- oxalate_olsen[,c(14,6)]
olsen_AL <- olsen_AL[,c(6,5)]

# run conversion function
converted_AL <- convert_P(nutr_df$pOxalate, oxalate_olsen, olsen_AL)

# add to dataframe
nutr_df$P_mean <- converted_AL$mean
nutr_df$P_sd <- converted_AL$sd


#### Convert potassium ####

# while we're waiting on the full dataset, here are averages we can use
# from Zebec et al: https://link.springer.com/content/pdf/10.1134%2FS1064229317130051.pdf
# these are in mg K/kg; divide by 390 to get them to mg/Kg 
AL <- c(179.57, 202.89, 210.26, 188.07, 197.05, 86.66, 197.98, 195.13, 184.98, 165.99)
AA <- c(170, 183.09, 197.23, 177.28, 168.88, 83.14, 187.73, 182.01, 190.9, 166.94)
AA_AL <- cbind(AA, AL)

# these are in mg K/kg; divide by 390 to get them to centimole/Kg (NEON units) 
# (according to https://www.dpi.nsw.gov.au/__data/assets/pdf_file/0007/127276/Chemical-tests.pdf)
AA_AL <- AA_AL/390

# run function to convert NEON's ammonium acetate values to Tedersoo's ammonium lactate
converted_K <- convert_K(nutr_df$kNh4d, AA_AL)

# add to dataframe
nutr_df$K_mean <- converted_K$mean
nutr_df$K_sd <- converted_K$sd


#### Get log-transformed weighted averages ####
# taking logs leads to negative values because our units are off (most likely). ignore for now.
# take the natural log of the nutrients
# nutr_df[,c("LogCa","LogK","LogMg","LogP")] <- log(nutr_df[,c("caNh4d","K_mean","mgNh4d","P_mean")]) 

# get means of log-transformed nutrient values, weighted by bulk density.
# nutr_tbl <- t(sapply(split(nutr_df, nutr_df$plotID.x), function(x) apply(x[, c("LogCa","LogK","LogMg","LogP")], 2, weighted.mean, 
#                                                                         x$bd))) 

# plot-level weighted averages of untransformed values
nutr_tbl <- t(sapply(split(nutr_df, nutr_df$plotID.x), function(x) apply(x[, c("caNh4d","K_mean","mgNh4d","P_mean")], 2, weighted.mean, 
                                                                         x$bd))) 

# remove any completely empty rows
nutr_tbl <- nutr_tbl[-which(rowSums(is.na(nutr_tbl))==4), ]

# check results dataframe
head(nutr_tbl) 

# compare the plot-level nutrient data we generated with the plots present in our microbial abundance data.
map <- readRDS(obs.table_16S.path) # read in observations
obs_plots <- unique(map$plotID) # get unique observation plotIDs
plots <- unique(rownames(nutr_tbl)) # get unique nutrient plotIDs

# um. we are missing nutrient data for most of the plots that we have observations for.
setdiff(obs_plots, plots) # plotIDs only in our observations: 48
setdiff(plots, obs_plots) # plotIDs only in the nutrient dataset
intersect(obs_plots, plots) # plotIDs common to both: 27
# DAMN IT

