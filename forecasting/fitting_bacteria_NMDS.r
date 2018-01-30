#Analyzing NMDS scores of bacterial composition
rm(list=ls())
library(vegan)

#load and proportionally normalize OTU table
otu <- readRDS('/fs/data3/caverill/NEFI_microbial/map_otu/16S_otu_clean.rds')
pro.function <- function(otu){
  for(i in 1:ncol(otu)){
    otu[,i] <- otu[,i] / sum(otu[,i])
  }
  return(otu)
}
otu <- pro.function(otu)

#get NMDS scores
bray.otu <- vegdist(t(otu), method = "bray")
nmds.otu <- metaMDS(bray.otu)
#check that it preserves the rank order well with the stress plot. It does.
stressplot(nmds.otu)
nmds.otu <- plot(nmds.otu)$sites

