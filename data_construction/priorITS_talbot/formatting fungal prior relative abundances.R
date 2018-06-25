#Building bray curtis similarity and predictor distance matrices for MRM analysis of Talbot et al. 2014
#first clear R environment.
rm(list=ls())
library(data.table)
library(vegan)
library(tidyr)
source('Scripts/space_time_functions.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/fg_assign.r')
source('/home/caverill/NEFI_microbe/NEFI_functions/lineaus.r')

#output save_path
output.path <- '/fs/data3/caverill/NEFI_microbial/prior_data/tal_all_prior_data.rds'

#number of genera to keep (top 10 or 20 most abundant)
n.gen <- 20

#load mapping file.
map <- readRDS('/fs/data3/caverill/Microbial_Space_Time_data/talbot_2014.data/Talbot2014_mapping_with_spatial.rds')
map <- data.table(map)
#load otu file
otu <- read.csv('/fs/data3/caverill/Microbial_Space_Time_data/talbot_2014.data/DOB_Soil_OTU_Table_11Feb13_UNITETAX_with_taxonomy.csv',header = T)

#load ecto hydrophobic status from hobbie
em.trait <- data.table(read.csv('/fs/data3/caverill/NEFI_microbial/ecto_genus_traits_hobbie_Jan2018.csv'))

#remove otu ID column, peel off taxonomy.
otu <- otu[,-1]
tax <- otu[,600]

#This is the bray dissimilarity matrix used for publication.
#subset mapping and OTU to only include samples present in this matrix.
bray.dis <- read.csv("/fs/data3/caverill/Microbial_Space_Time_data/talbot_2014.data/Talbot2014_brayCurtis500Avg_All.csv", header = TRUE, row.names = 1, check.names = FALSE)
otu <- otu[,colnames(otu) %in% colnames(bray.dis)] 
map <- map[map$Mapping.ID %in% colnames(otu)]

#sort so they are in the same order
otu <- otu[,sort(colnames(otu))]
sort.key <- colnames(otu)
map <- map[match(sort.key,map$Mapping.ID),]

#get day of year
#test <- strftime(map$Date.Sampled2, format = "%j")
map$doy <- as.Date(map$Date.Sampled,format='%m/%d/%Y')
map$year <- lubridate::year(map$doy)
map$doy  <- lubridate::yday(map$doy)
map$time <- ifelse(map$year == 11, map$doy,
                   ifelse(map$year == 12, map$doy + 365,
                          ifelse(map$year == 13, map$doy + 730, NA)))


#normalize the otu table
pro.function <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}
otu <- pro.function(otu)
#make sure column sums are 1.
colSums(otu)

#Get ECM, AM and SAP subsets of OTU table.
#for some reason 134 otus in this table are of the motif "k_Fungi, p_Basidiomycota...", instead of "k__Fungi, p_Basidiomycota..."
tax <- as.character(tax)
split.tax <- linaeus(tax)
for(i in 1:7){
  split.tax[,i] <- trimws(split.tax[,i])      #trim leading trailing white space.
}
#conditionally trimming.
for(j in 1:nrow(split.tax)){
  #if line is the "k__" motif, take everything from the 4th character on.
  if(split.tax[j,1]=='k__Fungi'){
    for(i in 1:ncol(split.tax)){
      split.tax[j,i] <- substring(split.tax[j,i],4)
    }
  }
  #if line is the "k_" motif, take everything from the 3rd character on.
  if(split.tax[j,1]=='k_Fungi'){
    for(i in 1:ncol(split.tax)){
      split.tax[j,i] <- substring(split.tax[j,i],3)
    }
  }
}

tax <- cbind(tax, split.tax)
#remove "possible" assignments


#Go ahead and assign function using FunGuild with the fg_assign function
tax <- fg_assign(tax)

#assign hydrophobic/hydrophillic
tax <- data.table(tax)
tax[genus %in% em.trait[ hydrophobic == 1,]$genus, hydrophobic := 1]
tax[genus %in% em.trait[hydrophillic == 1,]$genus,hydrophillic := 1]
tax[is.na( hydrophobic),  hydrophobic := 0]
tax[is.na(hydrophillic), hydrophillic := 0]


#get ten most abundant genera in dataset.
genera <- unique(tax$genus)
test <- data.table(cbind(tax, otu))
seq.out <- list()
for(i in 1:length(genera)){
  z <- test[genus == genera[i],]
  out <- colSums(z[,18:ncol(z)])
  seq.out[[i]] <- out
}
seq.out <- do.call('rbind',seq.out)

#Tried checking which genera are observed in every sample.
#This is only true for "unidentified".
#test <- data.frame(cbind(genera,seq.out))
#test <- test[!(apply(test, 1, function(y) any(y == 0))),]


#Count genus level abundance, grab some number of most abundant genera
counts <- rowSums(seq.out)
k <- data.table(cbind(genera,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(genera %in% c('unidentified'))] #remove the most abundant genus, "unidentified".
#grab genera of interest.
of_interest <- k$genera[1:n.gen]

#1. Get relative abundances of the most abundant genera
gen.list <- list()
for(i in 1:length(of_interest)){
  z <- data.table(cbind(tax,otu))
  z <- z[genus %in% of_interest[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  gen.list[[i]] <- out
}
gen.list <- data.frame(t(do.call('rbind',gen.list)))
colnames(gen.list) <- of_interest
gen.list$Mapping.ID <- rownames(gen.list)

#2. Get relative abundances of hydrophillic and hydrophobic taxa.
z <- data.table(cbind(tax,otu))
phil <- z[hydrophillic == 1,]
phob <- z[hydrophobic  == 1,]
start <- ncol(tax) + 1
out.phil <- colSums(phil[,start:ncol(phil)])
out.phob <- colSums(phob[,start:ncol(phob)])
hydro.out <- data.frame(t(rbind(out.phil,out.phob)))
colnames(hydro.out) <- c('hydrophillic','hydrophobic')
hydro.out$Mapping.ID <- rownames(hydro.out)

#3. get relative abundances of functional groups (ECM, AM, SAP, WR)
function_groups <- c('Ectomycorrhizal','Arbuscular','Saprotroph')
fun.list <- list()
for(i in 1:length(function_groups)){
  z <- data.table(cbind(tax,otu))
  z <- z[grep(function_groups[i],z$guild),]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  fun.list[[i]] <- out
}
fun.list <- data.frame(t(do.call('rbind',fun.list)))
colnames(fun.list) <- function_groups
fun.list$Mapping.ID <- rownames(fun.list)


#Merge together data you want for analysis.
abundances <- merge(fun.list,gen.list)
abundances <- merge(abundances,hydro.out)

#grab columns from map actually of interest.
map <- map[,.(Mapping.ID,Site,Plot,Horizon,longitude.dd,latitude.dd,pH,Perc.Soil.Moisture,Perc.N,Perc.C,CNRatio,Date.Sampled,doy,epoch.date,NPP,MAT,MAP,MAT_CV,MAP_CV)]
map <- merge(map,abundances)

#rename somethings
setnames(map,c('Perc.Soil.Moisture','Perc.N','Perc.C','CNRatio','latitude.dd','longitude.dd'),c('moisture','pN','pC','cn','latitude','longitude'))

#save output
saveRDS(map,output.path)