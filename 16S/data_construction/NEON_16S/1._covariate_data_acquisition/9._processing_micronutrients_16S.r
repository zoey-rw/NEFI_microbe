# Get micronutrient (Ca, Mg, P, K) data for each NEON plot
# from initial characterization data - mostly 2016
# NEON data is processed to better match the methods of the priors
# Ca and Mg - same methods, only units change.
# K - method changes and units change.
# P - NEON used two methods that are changed to match those of the priors.
# P and K will each have a mean, uncertainty, and variance.

rm(list=ls())
source('paths.r')
library(dplyr)

# read in soil chemistry and soil physical data
site_output_chem <- readRDS(dp1.10008_dp.10047_output_16S.path)[[1]]
site_output_phys <- readRDS(dp1.10008_dp.10047_output_16S.path)[[2]]

# read in raw data comparing P-extraction methods
p_methods <- read.csv(paste0(pecan_gen_16S_dir, "phosphorus_methods.csv"))

# merge nutrient data with bulk density data
merged <- merge(site_output_phys, site_output_chem, 
                by = c("horizonID"))

# drop unnecessary columns (see NEON variables files for information on each)
nutr_df <- merged[,(colnames(merged) %in% c("uid.x", "siteID.x", "plotID.x", "collectDate.x", "horizonID",
                                            "horizonName.x", "biogeoTopDepth", "biogeoBottomDepth",
                                            "mgNh4d", "caNh4d", "kNh4d", "OlsenPExtractable", "MehlichIIITotP", 
                                            "bulkDensTopDepth", "bulkDensBottomDepth", "bulkDensVolume", "bulkDensWetWeight",
                                            "bulkDensDryWeight", "bulkDensOvenDry", "bulkDensFieldMoist"))]

# calculate bulk density for the rows where it is empty
# dividing the dry weight (in grams) by the volume (in cm3)
# for some reason, rows that don't have calculated bulk density DO have volume + mass, and vice versa.
nutr_df$bd <- NA
nutr_df$bd <- ifelse(!is.na(nutr_df$bulkDensOvenDry), 
                     nutr_df$bulkDensOvenDry, 
                     nutr_df$bulkDensDryWeight/nutr_df$bulkDensVolume)

#### Convert phosphorus ####

# model relationship between MehlichIII P and Ammonium lactate (AL) P
# these match the models from Erikkson 2009
fitM3 <- lm(P.AL_color ~ P.M3colorom., data=p_methods)  # linear regression model on full data
coefM3 <- summary(fitM3)[[4]][2,1]
intM3 <- summary(fitM3)[[4]][1,1]
seM3 <- summary(fitM3)[[6]]

# function for converting values (there is probably a simpler way to do this)
P_M3.to.P_AL <- function(P_M3) {
  P_AL <- coefM3 * P_M3 + intM3
  return(P_AL)
}

# model relationship between Olsen P and Ammonium lactate (AL) P
fitO <- lm(P.ALICP ~ Olsen.PICP, data=p_methods)  # linear regression model on full data
coefO <- summary(fitO)[[4]][2,1]
intO <- summary(fitO)[[4]][1,1]
seO <- summary(fitO)[[6]]

# function for converting values
P_O.to.P_AL <- function(P_O) {
  P_AL <- 3.2778 * P_O + 27.88
  return(P_AL)
}

# initialize empty P columns
nutr_df$p <- NA
nutr_df$p_method <- NA

# set P to MehlichIII, and add source
nutr_df$p <- nutr_df$MehlichIIITotP
nutr_df[!is.na(nutr_df$p),]$p_method <- "MehlichIII"

# those that are still NA, set to Olsen P and add source
nutr_df[is.na(nutr_df$p),]$p <- nutr_df[is.na(nutr_df$p),]$OlsenPExtractable
nutr_df[is.na(nutr_df$p_method),]$p_method <- "OlsenP"

# convert P values
nutr_df$new_p <- NA
nutr_df[nutr_df$p_method == "OlsenP",]$new_p <- P_O.to.P_AL(nutr_df[nutr_df$p_method == "OlsenP",]$p)
nutr_df[nutr_df$p_method == "MehlichIII",]$new_p <- P_M3.to.P_AL(nutr_df[nutr_df$p_method == "MehlichIII",]$p)

# add P standard error
nutr_df$P_se <- NA
nutr_df[nutr_df$p_method == "OlsenP",]$P_se <- seO
nutr_df[nutr_df$p_method == "MehlichIII",]$P_se <- seM3


#### Convert potassium ####

# according to Zebec et al, K extracted using AL (priors) has a .97 correlation coeff with K extracted using AA (NEON)
# they don't give us raw data, but they do tell us n = 165, and that AA extracts 5% less then AL
# https://link.springer.com/content/pdf/10.1134%2FS1064229317130051.pdf

# what do I do with this information? how do we use the correlation coeff as uncertainty? 
# using the above values, the standard error of the correlation is .019 but that's probably not what we want
nutr_df$K <- nutr_df$kNh4d * .95

# take the natural log of the nutrients
nutr_df[,c("LogCa","LogK","LogMg","LogP")] <- log(nutr_df[,c("caNh4d","K","mgNh4d","new_p")]) 

# take the horizons that do not go deeper than 30 inches.
nutr_df <- nutr_df[nutr_df$bulkDensBottomDepth <= 30,]
# or take the top two horizons:
# nutr_df <- tbl_df(nutr_df) %>% group_by(plotID.x) %>% top_n(-2, bulkDensTopDepth)

# get means of log-transformed nutrient values, weighted by bulk density.
nutr_tbl <- t(sapply(split(nutr_df, nutr_df$plotID.x), function(x) apply(x[, c("LogCa","LogK","LogMg","LogP")], 2, weighted.mean, 
                                                                       x$bd))) 
# check results dataframe
head(nutr_tbl) 


# alright, negative values, that's not good. our units must be different from Tedersoo's, but they don't state units in their methods.
# NEON nutrients are in centimolesPerKilogram, other than phosphorus, which is milligramsPerKilogram
# check against bahram data - these values are already log-transformed
# bahram <- readRDS(bahram_metadata.path)

# compare the plot-level nutrient data we generated with the plots present in our microbial abundance data.
map <- readRDS(obs.table_16S.path) # read in observations
obs_plots <- unique(map$plotID) # get unique observation plotIDs
plots <- unique(nutr_tbl$plotID) # get unique nutrient plotIDs

# um. we are missing nutrient data for most of the plots that we have observations for.
setdiff(obs_plots, plots) # plotIDs only in our observations
setdiff(plots, obs_plots) # plotIDs only in the nutrient dataset
intersect(obs_plots, plots) # plotIDs common to both
# DAMN IT

# # multiply concentration by bulk density
# nutr_df$Ca_by_bd <- nutr_df$LogCa * nutr_df$bd
# nutr_df$K_by_bd <- nutr_df$LogK * nutr_df$bd
# nutr_df$Mg_by_bd <- nutr_df$LogMg * nutr_df$bd
# nutr_df$P_by_bd <- nutr_df$LogP * nutr_df$bd
# # sum the nutrient stocks from the top horizons, then divide by the sum of the volumes.
# nutr_tbl <- aggregate(cbind(nutr_df$Ca_by_bd, nutr_df$K_by_bd, nutr_df$Mg_by_bd, nutr_df$P_by_bd), 
#             by=list(Category=nutr_df$plotID.x), FUN=sum)
# colnames(nutr_tbl) <- c("plotID", "Ca", "K", "Mg", "P")
# nutr_tbl$volume <-  aggregate(nutr_df$bulkDensVolume, 
#                     by=list(Category=nutr_df$plotID.x), FUN=sum)