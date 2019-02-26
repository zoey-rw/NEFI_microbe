
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

