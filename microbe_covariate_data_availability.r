#Check with site by year_month combinations from mapping file even have data.
#assumes you have already ran 'formatting_NEON_microbial/1.cleaning_raw_map-otu_files.r'
#clear R environment, load packages, source functions.
rm(list=ls())
library(nneo)
library(data.table)
source('NEFI_functions/check_neon_products.r')

map <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_map_clean.rds')
of_interest <- c('DP1.10078.001', #soil chemical
                 'DP1.10086.001', #soil physical
                 'DP1.10100.001', #soil stable isotopes
                 'DP1.10104.001', #soil microbial biomass
                 'DP1.10080.001', #soil inorganic nitrogen pools and transformations
                 'DP1.10033.001', #litterfall and fine woody debris sampling
                 'DP1.10031.001', #litter chemical
                 'DP1.10026.001', #plant foliar physical and chemical properties
                 'DP1.10014.001', #coarse down wood bulk density sampling
                 'DP1.00098.001', #coarse down wood log survey
                 'DP1.10066.001', #root sampling (megapit)
                 'DP1.10102.001'  #root chemical properties
                 )

out <- check_neon_products(map,of_interest)
colSums(out[,4:ncol(out)])
