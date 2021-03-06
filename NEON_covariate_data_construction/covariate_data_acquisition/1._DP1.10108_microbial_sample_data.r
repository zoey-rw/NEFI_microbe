# Getting NEON soil core microbial metadata

rm(list=ls())
source('paths.r')
library(neonUtilities)

output.path <- dp1.10108.001_output.path

# Grab all sites that have microbial sequence metadata.
# Downloads about 12 MB, from 2013 to 2019, into the workspace. Takes a couple minutes.
dat <- loadByProduct(dpID="DP1.10108.001", site="all", startdate="2013-01", enddate="2019-05",
                     package="expanded", check.size=FALSE)

# Rename dataframes
dna.data <- dat$mmg_soilDnaExtraction
pcr.data_ITS <- dat$mmg_soilPcrAmplification_ITS
pcr.data_16S <- dat$mmg_soilPcrAmplification_16S

# Remove unnecessary columns from DNA data
dna.merge <- dna.data[,!(colnames(dat$mmg_soilDnaExtraction) %in% colnames(pcr.data_16S))]
keep_cols <- c("dnaSampleID", "namedLocation", "siteID", "plotID", "collectDate","deprecatedVialID")
dna.merge[,keep_cols] <- dna.data[,keep_cols]

# Merge 16S DNA/PCR data
core.data_16S <- merge(pcr.data_16S, dna.merge, all.x=T)
core.data_16S <- core.data_16S[which(core.data_16S$sequenceAnalysisType != "metagenomics"),]
core.data_16S$dateID <- substr(core.data_16S$collectDate, 1, 7)

# Merge ITS DNA/PCR data
core.data_ITS <- merge(pcr.data_ITS, dna.merge, all=T)
core.data_ITS <- core.data_ITS[which(core.data_ITS$sequenceAnalysisType != "metagenomics"),]
core.data_ITS$dateID <- substr(core.data_ITS$collectDate, 1, 7)


# Remove columns from 16S data, merge with ITS
core.data_16S.merge <- core.data_16S[,!(colnames(core.data_16S) %in% colnames(core.data_ITS))]
keep_cols <- c("geneticSampleID", "namedLocation", "siteID", "plotID", "collectDate", "dateID")
core.data_16S.merge[,keep_cols] <- core.data_16S[,keep_cols]
core.data_merged <- merge(core.data_ITS, core.data_16S.merge, all=T)


# Save DP1.10108 output
output <- list(core.data_ITS, core.data_16S, core.data_merged)
names(output) <- c("core_ITS", "core_16S", "core_16S_ITS")
saveRDS(output, output.path)


# get nested list of sites and dates sampled.
# for ITS samples
sites <- unique(output$core_ITS$siteID)
sites <- sites[!is.na(sites)]
site_dates_ITS <- list()
for(i in 1:length(sites)){
  dates <- unique(output$core_ITS[output$core_ITS$siteID == sites[i],]$dateID)
  dates <- dates[!is.na(dates)]
  site_dates_ITS[[i]] <- dates
}
names(site_dates_ITS) <- sites

# for 16S samples
sites <- unique(output$core_16S$siteID)
sites <- sites[!is.na(sites)]
site_dates_16S <- list()
for(i in 1:length(sites)){
  dates <- unique(output$core_16S[output$core_16S$siteID == sites[i],]$dateID)
  dates <- dates[!is.na(dates)]
  site_dates_16S[[i]] <- dates
}
names(site_dates_16S) <- sites
sites_dates_output <- list(site_dates_ITS, site_dates_16S)
names(sites_dates_output) <- c("sites_dates_ITS", "site_dates_16S")
 
#save site_dates output.
saveRDS(sites_dates_output,site_dates.path)
