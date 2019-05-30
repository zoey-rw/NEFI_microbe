#Getting NEON qPCR data for ITS and 16S data construction.
#outputs table of fungal, bacterial and archaeal qPCR data.
rm(list=ls())
source('paths.r')
library(neonUtilities)

# Specify output path
output.path <- dp1.10109.001_output.path

# Grab all sites that have soil core data
# Downloads about 1.6 MB, from 2013 to 2015, into the workspace. Takes a few minutes.
# Not downloading recent years because their qPCR protocol changed (universal 16S primer combines archaea and bacteria)
dat <- loadByProduct(dpID="DP1.10109.001", site="all", startdate="2013-01", enddate="2015-01",
                     package="expanded", check.size=T)

qPCR <- dat$mga_soilGroupAbundances

# Each unique sample ID is reported multiple times for each PCR target (ITS or 16S, and 16S can be setup to target bacteria or archaea).
# Would be better to query raw, but lets reformat this some so it pairs at the core scale better downstream.
fun <- qPCR[qPCR$targetTaxonGroup == 'fungi',]
names(fun)[names(fun) == 'copyNumberStandardDeviation'] <- 'fun.copyNumberStandardDeviation'
names(fun)[names(fun) == 'meanCopyNumber'] <- 'fun.CopyNumber'

bac <- qPCR[qPCR$targetTaxonGroup == 'bacteria',]
names(bac)[names(bac) == 'meanCopyNumber'] <- 'bac.CopyNumber'
names(bac)[names(bac) == 'copyNumberStandardDeviation'] <- 'bac.copyNumberStandardDeviation'
bac.merge <- bac[,c("dnaSampleID", "namedLocation", "siteID", "plotID", "collectDate","geneticSampleID", 'bac.CopyNumber','bac.copyNumberStandardDeviation')]

arc <- qPCR[qPCR$targetTaxonGroup == 'archaea',]
names(arc)[names(arc) == 'meanCopyNumber'] <- 'arc.CopyNumber'
names(arc)[names(arc) == 'copyNumberStandardDeviation'] <- 'arc.copyNumberStandardDeviation'
arc.merge <- arc[,c("dnaSampleID", "namedLocation", "siteID", "plotID", "collectDate","geneticSampleID", 'arc.CopyNumber','arc.copyNumberStandardDeviation')]

bac_arc <- merge(bac.merge,arc.merge, all = T)
qPCR.output <- merge(fun,bac_arc, all = T)

#save DP1.10109.00 output.
saveRDS(output, output.path)
