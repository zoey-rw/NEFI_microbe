##### preparing Bahram prior - taxa abundance, using taxonomic table and metadata. #####

# outputs 20 most cosmopolitan genera
# functional groups are not yet assigned

##### clear environment, source packages, functions and paths #####
rm(list=ls())

# devtools::install_version("lubridate", "1.6.0") # something about Red Hat Enterprise Linux 6; currently need to use older version
library(lubridate)  
library(data.table)
library(rgdal)
library(raster)

#library(SRAdb)
source('paths.r')
source('NEFI_functions/worldclim2_grab.r')
source('NEFI_functions/arid_extract.r')

# number of genera to keep (top 10 or 20 most abundant)
n.gen <- 20
incl_southern_lat <- T

if (incl_southern_lat == T) {
  gen_output.path <- cosmo_output_16S_south_lat.path
  phyla_output.path <- phyla_output_16S_south_lat.path
} else {
  gen_output.path <- cosmo_output_16S.path
  phyla_output.path <- phyla_output_16S.path
}

##### load files ####
dir <- data.dir 
# SRA RunInfo table from: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP021922
SRA <- read.csv(paste0(dir,"16S/scc_gen/SraRunTable.csv"))

# load tedersoo mapping file
map <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))
map <- data.table(map)

# load SV table as otu file
otu <- readRDS(bahram_dada2_SV_table.path)

# load taxonomy
tax <- readRDS(bahram_dada2_tax_table.path)

# load metadata from Bahram and Tedersoo - sent to Colin 
if (incl_southern_lat == T) {
  metadata <- readRDS(bahram_metadata_south_lat.path) 
} else { 
  metadata <- readRDS(bahram_metadata.path)
}

##### taxonomic and functional assignment #####

# remove leading "k__" in taxonomy.
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

# subset otu table and tax table to only include observations in map file
metadata$Run <- as.character(metadata$Run)
otu <- otu[rownames(otu) %in% metadata$Run,]
metadata <- metadata[metadata$Run %in% rownames(otu),]
# order OTU table to match the mapping file
otu <- otu[order(rownames(otu), metadata$Run),]

# for column names to be lower case.
tax <- as.data.frame(tax)
colnames(tax) <- tolower(colnames(tax))

# remove taxa that do not assign to a kingdom from tax and otu table.
tax <- tax[tax$kingdom == 'Bacteria'|tax$kingdom == 'Archaea',] # only removes ~500 counts
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]



#determine cosmopolitan genera.----
#condition 1: present in greater than 50% of samples.
#condition 2: top 10% sequence abundance (sensu Delgado-Barquez et al. 2018, Science)
genera <- unique(tax$genus)
genera <- as.character(genera)
test <- data.table(cbind(tax, t(otu)))
seq.out <- list()
cosmo.out <- list()
for(i in 1:length(genera)){
  z <- test[genus == genera[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  cosmo <- length(out[out > 0]) / length(out)
  seq.out[[i]] <- out
  cosmo.out[[i]] <- cosmo
}
seq.out <- do.call('rbind',seq.out)
cosmo.out <- do.call('rbind',cosmo.out)
j <- data.table(cbind(genera,cosmo.out))
colnames(j)[2] <- 'cosmo'
j$cosmo <- as.numeric(j$cosmo)
j <- j[order(-cosmo),]
counts <- rowSums(seq.out)
k <- data.table(cbind(genera,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(genera %in% c('unidentified'))] #remove the genus "unidentified".
k <- k[genera!=""&!is.na(genera),] #remove NA and empty genera
j <- j[!(genera %in% c('unidentified'))] #remove the genus "unidentified".
j <- j[genera!=""&!is.na(genera),] #remove NA and empty genera
head(cbind(j,k), 20)

#7 genera present in > 50% of samples are also in the top 10%.
#cosmo_genera <- j[cosmo >= 0.50,]$genera

# 20 genera present in > 94% of samples.
#cosmo_genera <- j[cosmo >= 0.94,]$genera
cosmo_genera <- colnames(readRDS(cosmo_output_16S.path)$abundances)[2:21]

#Get seq abundances of cosmo genera.----
gen.list <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(cosmo_genera)){
  z <- k[genus == cosmo_genera[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  gen.list[[i]] <- out
}
gen.list <- data.frame(t(do.call('rbind',gen.list)))
colnames(gen.list) <- cosmo_genera
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(gen.list)
gen.list <- cbind(other,gen.list)
gen.list <- list(gen.list,seq_total)
names(gen.list) <- c('abundances','seq_total')
gen.list$rel.abundances <- gen.list$abundances / gen.list$seq_total


#save output.----
saveRDS(gen.list, gen_output.path)
























#Get seq abundances of top phyla.
phyla <- unique(tax$phylum)
phyla <- as.character(phyla)
test <- data.table(cbind(tax, t(otu)))
seq.out <- list()
phyla.out <- list()
for(i in 1:length(phyla)){
  z <- test[phylum == phyla[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  cosmo_phyla <- length(out[out > 0]) / length(out)
  seq.out[[i]] <- out
  phyla.out[[i]] <- cosmo_phyla
}
seq.out <- do.call('rbind',seq.out)
phyla.out <- do.call('rbind',phyla.out)
j <- data.table(cbind(phyla,phyla.out))
colnames(j)[2] <- 'presence'
j$presence <- as.numeric(j$presence)
j <- j[order(-presence),]
counts <- rowSums(seq.out)
k <- data.table(cbind(phyla,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(phyla %in% c('unidentified'))] #remove the phylum "unidentified".
k <- k[phyla!=""&!is.na(phyla),] #remove NA and empty phyla
j <- j[!(phyla %in% c('unidentified'))] #remove the phylum "unidentified".
j <- j[phyla!=""&!is.na(phyla),] #remove NA and empty phyla
head(cbind(j,k), 15)

# get 15 most cosmopolitan phyla.
cosmo_phyla <- j$phyla[1:15]
cosmo_phyla <- j[presence >= 0.87,]$phyla

#Get seq abundances of cosmo genera.----
phyla.list <- list()
k <- data.table(cbind(tax,t(otu)))
for(i in 1:length(cosmo_phyla)){
  z <- k[phylum == cosmo_phyla[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  phyla.list[[i]] <- out
}
phyla.list <- data.frame(t(do.call('rbind',phyla.list)))
colnames(phyla.list) <- cosmo_phyla
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(phyla.list)
phyla.list <- cbind(other,phyla.list)
phyla.list <- list(phyla.list,seq_total)
names(phyla.list) <- c('abundances','seq_total')
phyla.list$rel.abundances <- phyla.list$abundances / phyla.list$seq_total

#save output.----
saveRDS(phyla.list, phyla_output.path)
