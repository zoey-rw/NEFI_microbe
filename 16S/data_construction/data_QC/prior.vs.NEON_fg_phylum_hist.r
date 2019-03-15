
rm(list=ls())
source('paths.r')

##### histograms of functional group abundances at NEON vs prior #####
all.neon <- readRDS(NEON_fg_abundances_16S.path)
all.bahram <- readRDS(prior_fg_abundances_16S.path)

par(mfrow=c(2,2))
for (p in 1:12) {
neon <- all.neon[[p]] # subset by functional group
bahram <- all.bahram[[p]]
rel.bahram <- as.data.frame(bahram$rel.abundances)[,2] # subset to relative abundance df
rel.neon <- as.data.frame(neon$rel.abundances)[,2]
hist(rel.bahram, main = paste("bahram", colnames(bahram$rel.abundances)[2])) # plot histogram
mtext(print(paste("mean rel.abundance:", round(mean(rel.bahram),3)))) # add mean
hist(rel.neon, main = paste("neon", colnames(neon$rel.abundances)[2]))
mtext(print(paste("mean rel.abundance:", round(mean(rel.neon, na.rm = T),3))))
}
dev.off()

##### histograms of phyla abundances #####
neon_phylo <- readRDS(NEON_16S_phylo_groups_abundances.path)
bahram_phylo <- readRDS(bahram_16S_common_phylo_groups_list.path)

n_phyla <- as.data.frame(neon_phylo[[1]][[2]]) # subset to phyla relative abundances
b_phyla <- as.data.frame(bahram_phylo[[1]][[2]])
n_phyla <- n_phyla[,order(colnames(n_phyla))] #alphabetize
b_phyla <- b_phyla[,order(colnames(b_phyla))]

par(mfrow=c(2,2))
for (p in 1:21) {
  rel.bahram <- b_phyla[,c(p)] # subset to rel.abundance values by phylum
  rel.neon <- n_phyla[,c(p)]
  hist(rel.bahram, main = paste("bahram", colnames(b_phyla)[p])) # plot histogram
  mtext(print(paste("mean rel.abundance:", round(mean(rel.bahram),3)))) # add means
  hist(rel.neon, main = paste("neon", colnames(n_phyla)[p]))
  mtext(print(paste("mean rel.abundance:", round(mean(rel.neon, na.rm = T),3))))
}


sv <- readRDS(NEON_dada2_SV_table_rare.path)
old_sv <- readRDS(NEON_dada2_SV_table.path)
phyla <- as.data.frame(neon_phylo[[1]][[2]])
phyla[phyla$Acidobacteria > .2,]

tax <- readRDS(NEON_dada2_tax_table.path)
dim(tax)
colnames(tax)

#format taxonomy table.----
#remove leading "k__".
for(i in 1:ncol(tax)){
  tax[,i] <- substring(tax[,i],4)
}

#for column names to be lower case.
colnames(tax) <- tolower(colnames(tax))

#remove taxa that do not assign to fungi from tax and otu table.----
tax <- as.data.frame(tax)

sv=old_sv; tax=tax; groups=unique(tax$phylum); tax_level="phylum"
phylum.out <- common_group_quantification(old_sv,tax,unique(tax$phylum),'phylum')
phylum.out$abundances <- as.data.frame(phylum.out$abundances)
phylum.out$abundances$rowsums <- rowSums(phylum.out$abundances, na.rm = T)

phylum.out$rel.abundances <- as.data.frame(phylum.out$rel.abundances)
bad <- rownames(phylum.out$rel.abundances[phylum.out$rel.abundances$Acidobacteria > .5,])
phylum.out$abundances[rownames(phylum.out$abundances) %in% bad,]

phylum.new <- common_group_quantification(sv,tax,unique(tax$phylum),'phylum')
phylum.new$abundances <- as.data.frame(phylum.new$abundances)
phylum.new$rel.abundances <- as.data.frame(phylum.new$rel.abundances)
phylum.new$abundances$rowsums <- rowSums(phylum.new$abundances, na.rm = T)


bad <- rownames(phylum.new$rel.abundances[phylum.new$rel.abundances$Acidobacteria > .4,])

phylum.new$abundances[phylum.new$abundances$rowsums < 10000,]
phylum.out$abundances[rownames(phylum.out$abundances) %in% bad,]




dat <- readRDS(hierarch_filled_16S.path)
core_mu <- dat$core.core.mu
dim(core_mu)
head(core_mu)

raw.truth <- readRDS(NEON_all.phylo.levels_plot.site_obs_16S.path)
raw.truth <- raw.truth$phylum
head(raw.truth)

map <- readRDS(core_obs_16S.path)

truth <- as.data.frame(raw.truth$core.fit)
truth$deprecatedVialID <- rownames(truth)
truth1 <- merge(truth, map[,c("deprecatedVialID", "geneticSampleID", "siteID")], by = "deprecatedVialID")
truth1$geneticSampleID <- gsub('-GEN','',truth1$geneticSampleID)
NEON <- merge(truth1, core_mu, by.x="geneticSampleID",by.y="sampleID")
dim(NEON)
plot(NEON$pH, NEON$Acidobacteria)


truth2[which(truth2$pH == truth2$pH[1]),]$siteID.x
class(truth2)
truth2$pH == truth2$pH[1]

meta <- readRDS(bahram_metadata.path)
dim(meta)
phyla <- bahram_phylo$phylum
phyla <- as.data.frame(phyla$rel.abundances)
phyla$Run <- rownames(phyla)

phyla <- phyla[order(phyla$Run),]
meta <- meta[order(meta$Run),]
plot(meta$pH, phyla$Acidobacteria)

ITS <- readRDS(hierarch_filled.path)
ITS <- ITS$core.core.mu
ITS[ITS$siteID=="HARV",]$pH


NEON[NEON$siteID.x=="CPER",]$pH
