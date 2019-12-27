# read in and combine all Ramirez synthesis data (3 CSV files from google drive)
source("paths_fall2019.r")
library(data.table)

output.path <- ramirez_raw_mapping_and_abundance.path

# define function for merging to replac
magic.merge <- function(x,y){
  require(dplyr)
  full_join(x,y) %>%
    inner_join(y) %>%
    coalesce(full_join(x,y) %>%
               inner_join(x))
}

file.dir <- file.path(scc_gen_16S_dir, "prior_abundance_mapping/Ramirez/")

# assuming file2 is the right one to use because it has all the cols/rows of file3, and more than file1.
file2 <- fread(file.path(file.dir, "merged_spec_env2.csv"), stringsAsFactors=FALSE, encoding="Latin-1")


# read in metadata for study 24
env_24 <- fread(file.path(file.dir, "env_24.csv"), stringsAsFactors=FALSE, encoding="Latin-1")

# fill in empty study numbers using row names
file2[is.na(file2$Study_refno),]$Study_refno <- as.numeric(gsub("X", "", file2[is.na(file2$Study_refno),]$dataset))
# fill in empty sample numbers using row names
file2[is.na(file2$Sample_Name_Study_refno2),]$Sample_Name_Study_refno2 <- file2[is.na(file2$Sample_Name_Study_refno2),]$Row.names
# fill in dataset names using Study_refno
file2[is.na(file2$dataset),]$dataset <- paste0("X", file2[is.na(file2$dataset),]$Study_refno)

# download name-matched file
name.matched <- curl::curl_download("https://static-content.springer.com/esm/art%3A10.1038%2Fs41564-017-0062-x/MediaObjects/41564_2017_62_MOESM7_ESM.zip", tempfile(fileext = ".zip"))
name.matched <- unzip(unzip(name.matched))
name.matched <- data.table::fread(name.matched[[1]])
name.matched <- as.data.frame(name.matched)
# subset to metadata.----
drop <- colnames(name.matched)[grepl("eukaryota|archaea|bacteria|unassigned", colnames(name.matched))]
metadata <- name.matched[,!colnames(name.matched) %in% drop]
# subset to taxonomy ----
bac <- colnames(name.matched)[grepl("p__bacteria", colnames(name.matched))]
bac <- name.matched[,colnames(name.matched) %in% bac]
colnames(bac) <- sapply(strsplit(colnames(bac), "__"), `[`, 3)
bac[, c("dataset","latitude","longitude","pH")] <- name.matched[,c("dataset","latitude","longitude","ph")]

# fix taxonomic names
file2_bac_names <- colnames(file2)[grepl("p_bacteria", colnames(file2))]
#not sure why this isn't working...
#file2_bac <- file2[which(colnames(file2) %in% file2_bac_names)]
file2_bac  <- subset(file2, select = colnames(file2) %in% file2_bac_names)
colnames(file2_bac) <- sapply(strsplit(colnames(file2_bac), "_"), `[`, 3)
file2_bac[, c("dataset","Sample_Name_Study_refno2")] <- file2[,c("dataset","Sample_Name_Study_refno2")]

# now let's use the taxonomy to merge file2 and the metadata from the taxonomy file. then use that metadata to get values from env24 file.
bac_merge <- bac[,c("acidobacteria", "actinobacteria", "bacteroidetes", "verrucomicrobia","dataset","latitude","longitude","pH")]
file2_bac_merge <- file2_bac[,c("acidobacteria", "actinobacteria", "bacteroidetes", "verrucomicrobia","dataset","Sample_Name_Study_refno2")]
tax_merged <- merge(file2_bac_merge, bac_merge)
tax_merged[tax_merged$dataset=="X24",]


colnames(tax_merged)[1:4] <- c("p_bacteria_acidobacteria", "p_bacteria_actinobacteria", "p_bacteria_bacteroidetes", "p_bacteria_verrucomicrobia")
file2$latitude <- as.numeric(file2$latitude)
file2$longitude <- as.numeric(file2$longitude)

df <- merge(file2, tax_merged, all.x=T, by = intersect(colnames(file2), colnames(tax_merged)))

# # fix a bunch of file type issues
# df$latitude <- as.numeric(df$latitude)
# df$longitude <- as.numeric(df$longitude)
# df$moisture <- as.numeric(df$moisture)
# df$pH <- as.numeric(df$pH)
# df$N <- as.numeric(df$N)
# df$C <- as.numeric(df$C)
# df$Study_refno <- as.numeric(df$Study_refno )

# 
# meta_merge <- metadata[,colnames(metadata) %in% c(intersect(colnames(metadata), colnames(df)))]
# df2 <- merge(df, meta_merge, by = intersect(colnames(meta_merge), colnames(df)), all.x = T)
# dim(df2)
# # df2[df2$dataset=="X24",1:80]
# df_with_24 <- merge(env_24[,c("Study_refno", "pH", "latitude","longitude", "C", "N", "moisture", "habitat","description","Natural_or_Managed")], df, by = c("Study_refno"), all.y=T)
# df_with_24[df_with_24$dataset=="X24",1:80]
# df_with_24[df_with_24$dataset=="X24"]$pH

# now merge other covariate data for study 24 into main df
env_24$dataset <- "X24"
env_merge <- env_24[,c("Study_refno", "pH", "latitude","longitude", "C", "N", "moisture", "habitat","description","Natural_or_Managed","sequening_platform")]
df_24 <- df[df$Study_refno=="24",]
out_24 <- magic.merge(env_merge, df_24)

df_no24 <- df[df$Study_refno!="24",]
out <- rbind(out_24, df_no24)


#meta_24 <- merge(env_24[,c("Study_refno", "pH", "latitude","longitude", "C", "N", "moisture", "habitat","description","Natural_or_Managed","sequening_platform")], df, by = c("latitude","longitude","pH"))

# merge back into main dataset
#df3 <- merge(file2, meta_24, by = c("Sample_Name_Study_refno2"), all.x=T)
# 
# df3[df3$dataset.x=="X24",1:80]
# out <- df3
# out[out$dataset.x=="X24",19600:19642]
# out$pH <- ifelse(!is.na(out$pH.x), out$pH.x, out$pH.y)
# out$dataset <- ifelse(!is.na(out$dataset.x), out$dataset.x, out$dataset.y)
# out$latitude <- ifelse(!is.na(out$latitude.x), out$latitude.x, out$latitude.y)
# out$longitude <- ifelse(!is.na(out$longitude.x), out$longitude.x, out$longitude.y)
# out$C <- ifelse(!is.na(out$C.x), out$C.x, out$C.y)
# out$N <- ifelse(!is.na(out$N.x), out$N.x, out$N.y)
# out$moisture <- ifelse(!is.na(out$moisture.x), out$moisture.x, out$moisture.y)
# out$habitat <- ifelse(!is.na(out$habitat.x), out$habitat.x, out$habitat.y)
# out$sequening_platform <- ifelse(!is.na(out$sequening_platform.x), out$sequening_platform.x, out$sequening_platform.y)


# get latitude and longitude back into this df
wanted_vars <- out[,c("Study_refno","Sample_Name_Study_refno2","dataset", "pH", "latitude","longitude", "C", "N", "moisture", "habitat")]
dim(wanted_vars)
dim(wanted_vars[complete.cases(wanted_vars),])
wanted_vars[!complete.cases(wanted_vars),]






forest_dats <- out[which(out$habitat %in% c("forest","Forest","logged forest")),]$dataset

# create list of studies to drop.----
# drop mouse poop samples
drop <- unique(out[which(out$habitat=="murine_stool"),]$dataset) #X6
# drop arctic and desert samples
drop <- append(drop, unique(out[which(out$biome.latitude.=="arctic"|out$biome.latitude.=="desert"),]$Study_refno))
# dropping X44 because it's data are wack (acknowledged in Ramirez et al.)
drop <- append(drop, "X44")
# XNA appears to be NEON....
drop <- append(drop, "XNA")

# add forest data.----
out$forest <- 0
out[which(out$dataset %in% forest_dats),]$forest <- 1

# drop any samples from our drop list.----
out <- out[!out$dataset %in% drop,]
dim(out)




saveRDS(out, output.path)