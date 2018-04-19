#making NEFI site predictions using JAGS models fit to Tedersoo data.
#clear environment, load packages.
rm(list=ls())
library(data.table)

#load JAGS model output.
ted.out <- readRDS('/fs/data3/caverill/NEFI_microbial/prior_data/ted_prior_JAGS_models.rds')

