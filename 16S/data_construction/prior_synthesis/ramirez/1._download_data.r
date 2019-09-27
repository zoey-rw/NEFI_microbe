## download relative abundance and site data from Ramirez et al. 16S synthesis ##
rm(list = ls())
library(data.table)
source('paths.r')
source('paths_fall2019.r')

# set output path
output.path <- ramirez_raw_mapping_and_abundance.path

# download name-matched file
name.matched <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM7_ESM.zip", tempfile(fileext = ".zip"))
name.matched <- unzip(unzip(name.matched))
name.matched <- data.table::fread(name.matched[[1]])
name.matched <- as.data.frame(name.matched)

# get dataset summary file
data.sum <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM5_ESM.csv", tempfile(fileext = ".csv"))
data.sum <- data.table::fread(data.sum)
data.sum <- as.data.frame(data.sum)


# subset to taxonomy metadata
drop <- colnames(name.matched)[grepl("eukaryota", colnames(name.matched))|grepl("archaea", colnames(name.matched))]
name.matched <- name.matched[,!colnames(name.matched) %in% drop]

# create study IDs
tax <- name.matched %>% group_by(dataset) %>% mutate(id = row_number())
tax <- as.data.frame(tax)
tax$sampleID <- paste(tax$dataset, tax$id, sep="-")
rownames(tax) <- tax$sampleID 
tax$id <- NULL

# remove unassigned and non-bacteria
metadata.cols <- colnames(tax)[!grepl("bacteria", colnames(tax))&!grepl("unassigned", colnames(tax))]
metadata <- tax[,colnames(tax) %in% metadata.cols]

tax <-  tax[,!colnames(tax) %in% metadata.cols]
tax$sampleID <- metadata$sampleID


# add biome etc from other data if possible
data.sum$dataset <- paste0("X", data.sum$OwnerRefNum.)
#to.match <- data.sum[,c("dataset","biome.latitude.","habitat","c","n","ph","latitude","longitude")]
# to.match <- data.sum[,c("biome.latitude.","habitat","latitude","longitude")]
# to.match <- to.match[!duplicated(to.match),]
# metadata <- data.frame(metadata)
# try1 <- merge(metadata, to.match, by=c("latitude","longitude"))
# dim(try1)
forest_dats <- data.sum[which(data.sum$habitat %in% c("forest","degraded forest")),]$dataset
# subset to unmanaged
#to.use <- mapping[which(mapping$natural_or_managed == "unmanaged" | mapping$natural_or_managed == "natural"),] 
#table(mapping$dataset)
to.use <- metadata

# create list of studies to drop:
# drop mouse poop samples
drop <- unique(data.sum[which(data.sum$habitat=="murine_stool"),]$dataset) #X6
# drop arctic and desert samples
drop <- append(drop, unique(data.sum[which(data.sum$biome.latitude.=="arctic"|data.sum$biome.latitude.=="desert"),]$dataset))
# dropping X44 because it's data are wack (acknowledged in Ramirez et al.)
drop <- append(drop, "X44")
# X1 for some reason, the values don't add up to anything near 1
drop <- append(drop, "X1")
# XNA appears to be NEON....
drop <- append(drop, "XNA")

# add forest data
metadata$forest <- 0
metadata[which(metadata$dataset %in% forest_dats),]$forest <- 1
to.use <- metadata[!metadata$dataset %in% drop,]
dim(to.use)
tax <- tax[which(tax$sampleID %in% to.use$sampleID),]

# remove "bacteria" column and add rownames
tax <- as.data.frame(tax)
tax$d__bacteria <- NULL
rownames(tax) <- tax$sampleID

output <- list(to.use, tax)

# save output.
saveRDS(output, output.path)
