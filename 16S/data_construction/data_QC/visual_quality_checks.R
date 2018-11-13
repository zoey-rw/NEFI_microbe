# creates NMDS plots of sample sequence depth and habitat.
# for interactive visual exploration. no saved outputs.
# not working yet - bray-curtis should be performed on the actual OTU table, 
# with sequence depth bin as a label. habitat can be a polygon or label to (hopefully) cluster samples.
# ideally samples WILL cluster in similarity by habitat, and not by sequence depth.

# load packages
library(vegan)

# read in sequence depth files
neon_track <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/NEON_processed/track.rds")
bahram_track <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/bahram_processed/track.rds")

#read in descriptors for samples & sites - just bahram for now
map<-readRDS("/fs/data3/caverill/NEFI_data/16S/pecan_gen/bahram_prior_gen.rds")

bahram_track <- as.data.frame(bahram_track)
bahram_track$Run <- rownames(bahram_track)

map_seqs <- merge(map, bahram_track[,c(4:5)], by="Run")
map_seqs$Seq_depth <- map_seqs$seqs_nonchim
map_seqs$Seq_bin <- as.numeric(cut(map_seqs$Seq_depth, 20))

Seq_depth <- map_seqs$Seq_depth
Seq_bin <- as.numeric(cut(map_seqs$Seq_depth, 20))
treat <- map_seqs$Biome

plot.df <- map_seqs[,.(Seq_depth, Seq_bin)]
dist.matrix<-as.matrix(plot.df)
rownames(dist.matrix) <- map_seqs$Run

################################
#### dissimilarity matrixes ####
################################


# try bray-curtis
vegdist(map_seqs$Seq_depth, method = "bray")
# try jaccard
vegdist(map_seqs$Seq_depth, method = "bray")

########################################################################
# plot NMDS
#######################################################################

all.MDS<-metaMDS(as.dist(log(dist.matrix +1)),k=2,zerodist="add")

#set up plot
all.site.coordinates<-data.frame(rownames(all.MDS))

#plot samples with Sample ID
ordiplot(all.MDS,type="text")
ordihull(all.MDS,groups=treat,draw="polygon",col="grey90",label=F)

#plot samples by seq depth
plot(all.MDS,type="n",display = "sites")
points(all.site.coordinates$NMDS1,all.site.coordinates$NMDS2,pch=c(21),bg=c("aliceblue","steelblue1","steelblue4","darkblue")[unclass(factor(Seq_bin))],cex=2.5)
ordiellipse(all.site.coordinates, Seq_bin, kind = "se", draw = "lines", conf = 0.95)
ordispider(all.MDS, Seq_bin, col = "blue", label= TRUE)



