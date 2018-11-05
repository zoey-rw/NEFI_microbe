{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf200
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;\red0\green0\blue0;\red27\green31\blue34;
\red0\green0\blue255;\red36\green38\blue41;\red104\green26\blue29;\red43\green39\blue19;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;\csgray\c0\c0;\cssrgb\c14118\c16078\c18039;
\cssrgb\c0\c0\c100000;\cssrgb\c18824\c20000\c21176;\cssrgb\c49020\c15294\c15294;\cssrgb\c22353\c20000\c9412;}
\margl1440\margr1440\vieww17760\viewh12540\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf2 \cb3 ##### preparing Bahram prior, using taxonomic table and metadata. #####\
\
# a lot of this is borrowed from the ITS script: \'934._format_tedersoo_2014_ITS_prior.r\'94\
# not completely functional yet - breaks before the most abundant genera are aggregated, \
# and functional groups are not yet assigned\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf2 ##### clear environment, source packages, functions and paths #####\
rm(list=ls())\
\
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf2 \expnd0\expndtw0\kerning0
# devtools::install_version("lubridate", "1.6.0") \kerning1\expnd0\expndtw0 # something about \expnd0\expndtw0\kerning0
Red Hat Enterprise Linux 6; \cf2 \cb3 \kerning1\expnd0\expndtw0 currently need to use older version\cf2 \cb3 \
library(lubridate)  \
library(data.table)\
\
#library(SRAdb)\
source('paths.r')\
source('NEFI_functions/worldclim2_grab.r')\
source('NEFI_functions/arid_extract.r')\
\
# number of genera to keep (top 10 or 20 most abundant)\
n.gen <- 20\
\
##### load files ####\
\
# SRA RunInfo table from: https://trace.ncbi.nlm.nih.gov/Traces/study/?acc=ERP021922\
SRA <- read.csv("/projectnb/talbot-lab-data/NEFI_data/16S/SraRunTable.csv")\
\
# metadata file that was sent to Colin\
metadata_raw <- read.csv('/projectnb/talbot-lab-data/NEFI_data/16S/metadata_bahram.csv')\
\
# from ITS script - ignore\
# map <- read.csv(ted_map_raw.path, header = TRUE, na.strings=c("", "NA"))\
# map <- data.table(map)\
\
# load SV table as otu file\
otu <- readRDS(bahram_dada2_SV_table.path)\
\
# load taxonomy\
tax <- readRDS(bahram_dada2_tax_table.path)\
\
# load times - sent to Colin separately by Leho Tedersoo\
time <- read.csv(ted_sampling_dates.path, header = TRUE, row.names=1, check.names = FALSE)\
\
\
##### Format metadata file #####\
\
# combine SRA sample names with Bahram sample names.\
sample_names <- SRA[,c(7,9,11)]\
\pard\pardeftab720\sl300\partightenfactor0
\cf2 \expnd0\expndtw0\kerning0
names(\cf2 \cb3 \kerning1\expnd0\expndtw0 sample_names\cf2 \cb3 \expnd0\expndtw0\kerning0
)[names(\cf2 \cb3 \kerning1\expnd0\expndtw0 sample_names\cf2 \cb3 \expnd0\expndtw0\kerning0
) == 'environment_biome'] <- 'Biome'\cf2 \cb3 \kerning1\expnd0\expndtw0 \
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0
\cf2 metadata <- merge(sample_names, metadata_raw, by.x = "Sample_Name", by.y = "X")\
metadata <- data.table(metadata)\cf2 \cb3 \
\
# subset to northern temperate latitudes\
metadata <- metadata[Lat < 66.5 & Lat > 23.5,]\
\
# format the time data frame (get rid of an empty column, etc.)\
colnames(time)[1] <- 'human.date'\
time[2] <- NULL\
#convert human readable date to days since epoch\
time$epoch.date <- round(\
  as.numeric(\
    as.POSIXlt(time$human.date, \
               format = "%m/%d/%y", origin = "01/01/1970"))/86400)\
# get day of year (doy) as well.\
time$doy  <- lubridate::yday(as.Date(time$human.date,format='%m/%d/%Y'))\
time$year <- lubridate::year(as.Date(time$human.date,format='%m/%d/%Y'))\
time$epoch.date <- (time$year - 9)*365 + time$doy\
\
# drop samples in time table not present in mapping file.\
time <- time[row.names(time) %in% metadata$Sample_Name,]\
\
# push times into the mapping file\
time$Sample_Name <- rownames(time)\
metadata <- merge(metadata,time, by = 'Sample_Name', all.x=T)\
\
# assign forest variables\
# shouldn't this also include  'Moist tropical forests', 'Tropical montane forests' and 'Southern_temperate_forests'?\
metadata$forest <-ifelse(metadata$Biome %in% c('Temperate_coniferous_forests','Temperate_deciduous_forests','Dry_tropical_forests','Boreal_forests'),1,0)\
metadata$conifer <- ifelse(metadata$Biome %in% c('Temperate_coniferous_forests'),1,0)\
\
##### taxonomic and functional assignment #####\
\
# remove leading "k__" in taxonomy.\
for(i in 1:ncol(tax))\{\
  tax[,i] <- substring(tax[,i],4)\
\}\
\
# subset otu table and tax table to only include observations in map file\
metadata$Run <- as.character(metadata$Run)\
otu <- otu[rownames(otu) %in% metadata$Run,]\
metadata <- metadata[metadata$Run %in% rownames(otu),]\
# order OTU table to match the mapping file\
otu <- otu[order(rownames(otu), metadata$Run),]\
\
# for column names to be lower case.\
tax <- as.data.frame(tax)\
colnames(tax) <- tolower(colnames(tax))\
\
# remove taxa that do not assign to bacteria or archaea from tax and otu table.\
tax <- tax[!is.na(tax$kingdom),]\
otu <- otu[,colnames(otu) %in% rownames(tax)]\
tax <- tax[rownames(tax) %in% colnames(otu),]\
\
# normalize the otu table \
otu <- t(otu)\
pro.function <- function(otu)\{\
  for(i in 1:ncol(otu))\{\
    otu[,i] <- otu[,i] / sum(otu[,i])\
  \}\
  return(otu)\
\}\
otu <- pro.function(otu)\
\
# make sure column sums are 1\
colSums(otu)\
\
# aggregate important classes and genera\
# get most abundant genera in dataset.\
genera <- unique(tax$genus)\
test <- data.table(cbind(tax, otu))\
seq.out <- list()\
for(i in 1:length(genera))\{\
  z <- test[genus == genera[i],]\
  start <- ncol(tax) + 1\
  out <- colSums(z[,start:ncol(z)])\
  seq.out[[i]] <- out\
\}\
\
# something breaks here\
# Count genus level abundance, grab some number of most abundant genera\
seq.out <- do.call('rbind',seq.out)\
counts <- rowSums(seq.out)\
k <- data.table(cbind(genera,counts))\
k$counts <- as.numeric(as.character(k$counts))\
k <- k[order(-counts),]\
k <- k[!(genera %in% c('unidentified'))] #remove the genus "unidentified".\
#grab genera of interest.\
of_interest <- k$genera[1:n.gen]\
\
# Get relative abundances of the most abundant genera\
gen.list <- list()\
for(i in 1:length(of_interest))\{\
  z <- data.table(cbind(tax,otu))\
  z <- z[genus %in% of_interest[i],]\
  start <- ncol(tax) + 1\
  out <- colSums(z[,start:ncol(z)])\
  gen.list[[i]] <- out\
\}\
gen.list <- data.frame(t(do.call('rbind',gen.list)))\
colnames(gen.list) <- of_interest\
gen.list$Mapping.ID <- rownames(gen.list)\
\
\
# Merge together data aggregated groups you want for analysis.\
# abundances <- merge(fun.list,gen.list)\
# abundances <- gen.list\
\
# final merging of files and worldclim grab.----\
# grab columns from map actually of interest.\
# metadata <- metadata[,.(tedersoo.code,SRR.id,Site,longitude,latitude,pH,Moisture,N,C,C_N,human.date,doy,epoch.date,NPP,forest,conifer,relEM)]\
# metadata <- merge(metadata,abundances, by.x = 'Run',by.y = 'Mapping.ID')\
# rename columns\
# setnames(metadata,c('tedersoo.code','Moisture','N' ,'C' ,'C_N'),\
         c('Mapping.ID','moisture','pN','pC','cn'))\
\
# save output.----\
# saveRDS(metadata, "/projectnb/talbot-lab-data/NEFI_data/16S/scc_gen/")\
}