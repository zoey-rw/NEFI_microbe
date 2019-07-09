# Download all ITS marker gene sequences from 5 NEON sites
# does not download early/provisional data - that was obtained directly from Lee Stanish
rm(list=ls())
library(RCurl)
library(dada2)
library(neonUtilities)
library(utils)
source('paths.r')
source('NEFI_functions/merge_left.r')
source('NEFI_functions/tic_toc.r')

site.list <- c("DSNY", "STER", "CPER", "OSBS", "HARV")

# set sequence and output directory
seq.dir <- paste0(big_data_dir, 'NEON_2014-2017/ITS/')
download.dir <- paste0(seq.dir, 'raw_fastq_ITS/')

# First, load metadata to determine the download URLS
its_metadata.path <- "/projectnb/talbot-lab-data/NEFI_data/big_data/NEON_2014-2017/sequence_metadata.rds" # delete this line later
if (file.exists(its_metadata.path)) { # Check existence of metadata
  its_metadata <- readRDS(its_metadata.path)
} else { # Download if needed
  its_metadata <- its_metadata.path
  # Get metadata for all NEON sites (function isn't vectorized, so this is easiest for now)
  dat <-  loadByProduct(dpID = "DP1.10108.001", site = "all", startdate = "2013-01", enddate = "2019-05", package = "basic", check.size = TRUE)
  its_mgs <- dat$mmg_soilMarkerGeneSequencing_ITS # ITS sequencing files
  its_pcr <- dat$mmg_soilPcrAmplification_ITS # ITS PCR files
  its_dna <- dat$mmg_soilDnaExtraction[dat$mmg_soilDnaExtraction$sequenceAnalysisType == "marker gene", ] # DNA extraction files
  its_pcr_mgs <- merge_left(its_mgs, its_pcr) # merge together using custom function. 
  its_pcr_mgs_dna <- merge_left(its_pcr_mgs, its_pcr, all.x = T) # "all.x = T": keep sample even if no DNA data matches
  its_pcr_mgs_dna <- its_pcr_mgs_dna[!duplicated(its_pcr_mgs_dna$dnaSampleID),] # remove duplicates
  its_metadata <- its_pcr_mgs_dna[its_pcr_mgs_dna$siteID %in% site.list, ] # subset to our 5 sites
  saveRDS(its_metadata, its_metadata.path) # save!
}

# only download files that we don't have
all_seq_runIDs <- unique(its_metadata$sequencerRunID)
downloaded.files <- list.files(download.dir) # get list of tarfiles we already have downloaded
downloaded_seq_runIDs <- unique(substr(downloaded.files, 5, 9)) # subset just the runID
not_downloaded <- setdiff(all_seq_runIDs, downloaded_seq_runIDs) # which tarfiles are we missing?
# there are two IDs like "150225_M02149_0120_000000000-ACC8C" - those were for provisional data, ignore

# download missing files
for (u in 1:length(not_downloaded)){ # loop over every undownloaded sequencer run
  
  potential_urls <- list() # create every possible url that this file could live in
  years <- c("2017", "2018", "2019", "2020") 
  for (y in 1:length(years)) { # any processing year, and R1 + R2
    year <- years[y]
    potential_urls <- append(potential_urls, paste0("https://neon-microbial-raw-seq-files.s3.data.neonscience.org/", year, "/BMI_", not_downloaded[u], "_ITS_R1.fastq.tar.gz"))
    potential_urls <- append(potential_urls, paste0("https://neon-microbial-raw-seq-files.s3.data.neonscience.org/", year, "/BMI_", not_downloaded[u], "_ITS_R2.fastq.tar.gz"))
  }
  
  for(i in 1:length(potential_urls)) { # loop over all potential urls
    url <- as.character(potential_urls[i])
    if(url.exists(url)){ # if it exists, download it
      cat(paste0("Downloading file: ", url, " \nto directory: ", download.dir))
      tic() # start timer
      download.file(url=url, destfile = paste0(download.dir, basename(url)))
      toc() # report how long the download took
    } # close download-existing-urls loop
  } # close check-potential-url-download loop
} # close undownloaded-sequence loop

# not relevant for our 5 sites, but:
# "BMI_C3CN4_ITS.fastq.tar.gz" is incorrectly listed as the processedSeqFileName for 120 samples 
# the internal lab ID would have to be used to link them to files

