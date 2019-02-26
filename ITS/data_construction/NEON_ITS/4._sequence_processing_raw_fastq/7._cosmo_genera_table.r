#NEON taxonomic group tables.
#common generea and phyla level.
#Phyla: Basidiomycota, Ascomycota, Zygomycota, Chytridmycota
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')

#set output paths.----
output.path1 <- NEON_ITS_fastq_cosmo_genera.path
output.path2 <- NEON_ITS_fastq_phyla.path

#load data.----
sv <- readRDS(NEON_ITS_fastq_SV.table_clean.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_ITS_fastq_fun_clean.path)
cosmo_ref <- readRDS(tedersoo_ITS_cosmo_genera_list.path)
cosmo_ref <- colnames(cosmo_ref$abundances)
cosmo_ref <- cosmo_ref[!cosmo_ref %in% c('other')]
phyla_ref <- c('Basidiomycota','Ascomycota')

#check that things are properly ordered.----
if(sum(rownames(fun) == colnames(sv)) != ncol(sv)){
  cat('Stop. rownames of tax-fun table do not match column names of SV table.\n')
}

#aggregate cosmopolitan genera.----
test <- data.table(cbind(fun, t(sv)))
seq.out <- list()
cosmo.out <- list()
for(i in 1:length(cosmo_ref)){
  z <- test[genus == cosmo_ref[i],]
  start <- ncol(fun) + 1
  out <- colSums(z[,start:ncol(z)])
  cosmo <- length(out[out > 0]) / length(out)
  seq.out[[i]] <- out
  cosmo.out[[i]] <- cosmo
}
genera.seq.out <- do.call('cbind',seq.out)
cosmo.out <- do.call('rbind',cosmo.out)
colnames(genera.seq.out) <- cosmo_ref

#aggregate common phyla.----
test <- data.table(cbind(fun, t(sv)))
phyla_seq.out <- list()
phyla.out <- list()
for(i in 1:length(phyla_ref)){
  z <- test[phylum == phyla_ref[i],]
  start <- ncol(fun) + 1
  out <- colSums(z[,start:ncol(z)])
  phylum <- length(out[out > 0]) / length(out)
  phyla_seq.out[[i]] <- out
  phyla.out[[i]] <- phylum
}
phyla_seq.out <- do.call(cbind,phyla_seq.out)
colnames(phyla_seq.out) <- phyla_ref
phyla.out <- do.call(rbind,phyla.out)


#get a handle on most abundant and most cosmopolitan genera.----
j <- data.table(cbind(cosmo_ref,cosmo.out))
colnames(j)[2] <- 'cosmo'
j$cosmo <- as.numeric(j$cosmo)
j <- j[order(-cosmo),]
counts <- colSums(genera.seq.out)
k <- data.table(cbind(cosmo_ref,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(cosmo %in% c('unidentified'))] #remove the genus "unidentified".
j <- j[!(cosmo %in% c('unidentified'))] #remove the genus "unidentified".
#head(cbind(j,k), 30)


#Get together genera abundances and seq_total.----
seq_total <- rowSums(sv)
seq_total <- seq_total[order(match(names(seq_total),rownames(genera.seq.out)))]
other <- seq_total - rowSums(genera.seq.out)
genera.abundances <- data.frame(other,genera.seq.out)
genera.rel.abundances <- genera.abundances / seq_total

for(i in 1:ncol(genera.abundances)){
  dat <- genera.abundances[,i]
  report <- paste0(round((length(dat[dat > 0]) / length(dat))*100,2), '% of samples have ',colnames(genera.abundances)[i],'\n')
  cat(report)
}

#Get together phyla abunadnces and seq_total.----
other <- seq_total - rowSums(phyla_seq.out)
phyla.abundances <- data.frame(other, phyla_seq.out)
phyla.rel.abundances <- phyla.abundances / seq_total

#get other IDs in here. deprecatedVialID does not match all products.----
    genera.abundances$geneticSampleID <- rownames(    genera.abundances)
genera.rel.abundances$geneticSampleID <- rownames(genera.rel.abundances)
     phyla.abundances$geneticSampleID <- rownames(     phyla.abundances)
 phyla.rel.abundances$geneticSampleID <- rownames( phyla.rel.abundances)

#save output.----
#cosmo genera.
genera.out <- list(genera.abundances,genera.rel.abundances,seq_total)
 phyla.out <- list( phyla.abundances, phyla.rel.abundances,seq_total)
names(genera.out) <- c('abundances','rel.abundances','seq_total')
names( phyla.out) <- c('abundances','rel.abundances','seq_total')
saveRDS(genera.out, output.path1)
saveRDS( phyla.out, output.path2)