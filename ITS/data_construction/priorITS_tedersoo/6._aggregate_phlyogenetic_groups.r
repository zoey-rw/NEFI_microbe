#Aggregate cosmpolitan genera.
#clear environment, source paths, packages and functions.
rm(list=ls())
source('paths.r')
library(data.table)

#set output paths.----
cosmo_output.path <- tedersoo_ITS_cosmo_genera_list.path
phyla_output.path <- tedersoo_ITS_phyla_list.path

#load data.----
map <- readRDS(tedersoo_ITS_clean_map.path)
otu <- readRDS(ted_2014_SV.table.path)
tax <- readRDS(ted_2014_tax.path)

#subset otu file to match mapping file.----
otu <- otu[rownames(otu) %in% map$SRR.id,]

#format taxonomy table.----
#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#for column names to be lower case.
colnames(tax) <- tolower(colnames(tax))

#remove taxa that do not assign to fungi from tax and otu table.
tax <- tax[tax$kingdom == 'Fungi',]
otu <- otu[,colnames(otu) %in% rownames(tax)]
tax <- tax[rownames(tax) %in% colnames(otu),]

#determine cosmopolitan genera.----
#condition 1: present in greater than 50% of samples.
#condition 2: top 10% sequence abundance (sensu Delgado-Barquez et al. 2018, Science)
genera <- unique(tax$genus)
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
j <- j[!(genera %in% c('unidentified'))] #remove the genus "unidentified".
head(cbind(j,k), 20)

#7 genera present in > 50% of samples are also in the top 10%.
cosmo_genera <- j[cosmo >= 0.50,]$genera

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
gen.list$rel.abundaces <- gen.list$abundances / gen.list$seq_total


#get abundances of phyla (Ascomycota, Basidiomycota).----
phyla.list <- list()
k <- data.table(cbind(tax,t(otu)))
phyla_ref <- c('Ascomycota','Basidiomycota')
for(i in 1:length(phyla_ref)){
  z <- k[phylum == phyla_ref[i],]
  start <- ncol(tax) + 1
  out <- colSums(z[,start:ncol(z)])
  phyla.list[[i]] <- out
}
phyla.list <- data.frame(t(do.call('rbind',phyla.list)))
colnames(phyla.list) <- phyla_ref
seq_total <- colSums(k[,start:ncol(k)])
other <- seq_total - rowSums(phyla.list)
phyla.list <- cbind(other,phyla.list)
phyla.list <- list(phyla.list,seq_total)
names(phyla.list) <- c('abundances','seq_total')
phyla.list$rel.abundaces <- phyla.list$abundances / phyla.list$seq_total


#save output.----
saveRDS(gen.list, cosmo_output.path)
saveRDS(phyla.list, phyla_output.path)
cat('YOU DID IT. GREAT JOB. <3 Colin.')
