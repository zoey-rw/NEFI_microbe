#Need link everything to IDs in covariate dataset, drop singletons, filter out samples w/ less than 1000 reads.
#As we drop samples we also need to update the tax/fun products as well.
#clear environment, source paths, packages and functions.
rm(list=ls())
library(data.table)
source('paths.r')

#set output paths.----
sv_output.path <- NEON_ITS_fastq_SV.table_clean.path
fun_output.path <- NEON_ITS_fastq_fun_clean.path

#load data.----
sv <- readRDS(NEON_ITS_fastq_SV.table.path)
map <- readRDS(hierarch_filled.path)
map <- map$core.obs
fun <- readRDS(NEON_ITS_fastq_fun.path)
tax <- readRDS(NEON_ITS_fastq_tax.path)

#subset to sv rows in map file. deprecatedVialID links map to sv.table.----
rownames(sv) <- toupper(rownames(sv))
sv.id <- rownames(sv)
map <- map[map$geneticSampleID %in% sv.id,]
sv <-  sv[rownames(sv) %in% map$geneticSampleID,]
rownames(map) <- map$geneticSampleID

#put in same order
map <- map[order(rownames(map)),]
sv <-  sv[order(rownames( sv)),]

#remove samples with less than 200 reads from map and sv.----
map$seq.depth <- rowSums(sv)
map <- map[map$seq.depth > 200,]
sv <- sv[rownames(sv) %in% rownames(map),]

#kill SVs that no longer have any sequences or are now singletons.----
sv.filter <- data.frame(colnames(sv),colSums(sv))
colnames(sv.filter) <- c('sv.ID','n.seqs')
to_remove <- sv.filter[sv.filter$n.seqs < 2,]$sv.ID
sv <-  sv[,!(colnames(sv) %in% to_remove) ]
fun <- fun[!(rownames(fun) %in% to_remove),]

#remove taxa that do not assign to fungi.----
to_remove <- rownames(fun[is.na(fun$kingdom),])
fun <- fun[!(rownames(fun) %in% to_remove),]
sv <-  sv[,!(colnames(sv) %in% to_remove) ]

#save output.----
saveRDS( sv,  sv_output.path)
saveRDS(fun, fun_output.path)
