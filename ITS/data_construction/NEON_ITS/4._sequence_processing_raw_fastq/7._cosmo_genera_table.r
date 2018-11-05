#NEON functional group tables.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')

#set output paths.----
output.path <- NEON_ITS_fastq_cosmo_genera.path

#load data.----
sv <- readRDS(NEON_ITS_fastq_SV.table_clean.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_ITS_fastq_fun_clean.path)
cosmo_ref <- readRDS(tedersoo_ITS_cosmo_genera_list.path)
cosmo_ref <- colnames(cosmo_ref$abundances)
cosmo_ref <- cosmo_ref[!cosmo_ref %in% c('other')]

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
seq.out <- do.call('cbind',seq.out)
cosmo.out <- do.call('rbind',cosmo.out)
colnames(seq.out) <- cosmo_ref

#get a handle on most abundant and most cosmopolitan genera.----
j <- data.table(cbind(cosmo_ref,cosmo.out))
colnames(j)[2] <- 'cosmo'
j$cosmo <- as.numeric(j$cosmo)
j <- j[order(-cosmo),]
counts <- colSums(seq.out)
k <- data.table(cbind(cosmo_ref,counts))
k$counts <- as.numeric(as.character(k$counts))
k <- k[order(-counts),]
k <- k[!(cosmo %in% c('unidentified'))] #remove the genus "unidentified".
j <- j[!(cosmo %in% c('unidentified'))] #remove the genus "unidentified".
#head(cbind(j,k), 30)


#Get together abundances and seq_total.----
seq_total <- rowSums(sv)
seq_total <- seq_total[order(match(names(seq_total),rownames(seq.out)))]
other <- seq_total - rowSums(seq.out)
abundances <- data.frame(other,seq.out)
rel.abundances <- abundances / seq_total

for(i in 1:ncol(abundances)){
  dat <- abundances[,i]
  report <- paste0(round((length(dat[dat > 0]) / length(dat))*100,2), '% of samples have ',colnames(abundances)[i],'\n')
  cat(report)
}

#get other IDs in here. deprecatedVialID does not match all products.----
abundances$deprecatedVialID <- rownames(    abundances)
rel.abundances$deprecatedVialID <- rownames(rel.abundances)

abundances <-     merge(    abundances, map[,c('deprecatedVialID','geneticSampleID')], 
                            by.x = 'deprecatedVialID', by.y = 'deprecatedVialID', all.x = T)
rel.abundances <- merge(rel.abundances, map[,c('deprecatedVialID','geneticSampleID')], 
                        by.x = 'deprecatedVialID', by.y = 'deprecatedVialID', all.x = T)

#save output.----
dat.out <- list(abundances,rel.abundances,seq_total)
names(dat.out) <- c('abundances','rel.abundances','seq_total')
saveRDS(dat.out, output.path)
