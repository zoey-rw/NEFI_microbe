# combine Ramirez synthesis data from google drive and paper supplement.
source("paths_fall2019.r")
library(data.table)
file.dir <- file.path(scc_gen_16S_dir, "prior_abundance_mapping/Ramirez/")

# first, download name-matched file of taxon abundances
name.matched <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM7_ESM.zip", tempfile(fileext = ".zip"))
name.matched <- unzip(unzip(name.matched))
name.matched <- data.table::fread(name.matched[[1]])
name.matched <- as.data.frame(name.matched)

# now read in file2, from which we will grab the sample names
file2 <- fread(file.path(file.dir, "merged_spec_env2.csv"), stringsAsFactors=FALSE, encoding="Latin-1")

# fill in empty study numbers using row names
file2[is.na(file2$Study_refno),]$Study_refno <- as.integer(gsub("X", "", file2[is.na(file2$Study_refno),]$dataset))
# fill in empty sample numbers using row names
file2[is.na(file2$Sample_Name_Study_refno2),]$Sample_Name_Study_refno2 <- file2[is.na(file2$Sample_Name_Study_refno2),]$Row.names
# fill in dataset names using Study_refno
file2[is.na(file2$dataset),]$dataset <- paste0("X", file2[is.na(file2$dataset),]$Study_refno)

file2_merge <- file2[,c("p_bacteria_planctomycetes","p_bacteria_proteobacteria","p_bacteria_chloroflexi","p_bacteria_actinobacteria","p_bacteria_armatimonadetes","p_bacteria_bacteroidetes","dataset","Sample_Name_Study_refno2","pH","habitat","latitude","longitude","biome.latitude.","C","N")]
file2_merge <- file2_merge[!duplicated(file2_merge[,1:7])]

# use taxon abundances to get the sample names.
tax_merged <- merge(file2_merge, name.matched,by.x=c("p_bacteria_planctomycetes","p_bacteria_proteobacteria","p_bacteria_chloroflexi","p_bacteria_actinobacteria","p_bacteria_armatimonadetes","p_bacteria_bacteroidetes","dataset"), by.y=c("p__bacteria__planctomycetes","p__bacteria__proteobacteria","p__bacteria__chloroflexi","p__bacteria__actinobacteria","p__bacteria__armatimonadetes","p__bacteria__bacteroidetes","dataset"))

# fix things that couldn't merge due to NA's.
df <- tax_merged
df$pH <- ifelse(is.na(df$ph), df$pH, df$ph)
df$latitude <- as.numeric(ifelse(is.na(df$latitude.y), df$latitude.x, df$latitude.y))
df$longitude <- as.numeric(ifelse(is.na(df$longitude.y), df$longitude.x, df$longitude.y))
df$sampleID <- df$Sample_Name_Study_refno2
#df$pC <- ifelse(is.na(df$longitude.y), df$longitude.x, df$longitude.y)

# now pull in the dataset summary file to get all the other predictors
data.sum <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM5_ESM.csv", tempfile(fileext = ".csv"))
data.sum <- data.table::fread(data.sum)
data.sum[1:5, 1:20] #check metadata - has sample names, dataset names, and full set of covariates (including carbon, nitrogen)
data.sum$dataset <- paste0("X", data.sum$study_refno)
master.df <- merge(df, data.sum, by=c("latitude","longitude","ph","dataset"))
master.df <- master.df[!duplicated(master.df[,c("latitude","longitude","ph","dataset")]),]
master.df$sampleID <- master.df$Sample_Name_Study_refno2
#forest_dats <- data.sum[which(data.sum$habitat %in% c("forest","Forest","logged forest")),]$dataset


out <- master.df
out$habitat <- ifelse(is.na(out$habitat.y), out$habitat.x, out$habitat.y)
forest_dats <- out[which(out$habitat %in% c("forest","Forest","logged forest")),]$dataset
# create list of studies to drop.----
# drop mouse poop samples
drop <- unique(out[which(out$habitat=="murine_stool"),]$dataset) #X6
# drop arctic and desert samples
drop <- append(drop, unique(out[which(out$biome.latitude..y=="arctic"|out$biome.latitude..y=="desert"),]$Study_refno))
# dropping X44 because it's data are wack (acknowledged in Ramirez et al.)
drop <- append(drop, "X44")
# XNA appears to be NEON....
drop <- append(drop, "XNA")

# add forest data.----
out$forest <- 0
out[which(out$dataset %in% forest_dats),]$forest <- 1
out[which(is.na(out$habitat)),]$forest <- NA

# drop any samples from our drop list.----
out <- out[!out$dataset %in% drop,]

# check that things worked alright.
plot(out$p__bacteria__acidobacteria, out$ph)
wanted_vars <- out[,c("dataset","Sample_Name_Study_refno2","dataset", "pH", "latitude","longitude", "c", "n", "moisture", "habitat","forest")]
dim(wanted_vars[complete.cases(wanted_vars),])
wanted_vars[!complete.cases(wanted_vars),]

saveRDS(out, ramirez_raw_mapping_and_abundance.path)
