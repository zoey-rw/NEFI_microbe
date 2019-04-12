# creates NMDS plots of sample sequence depth and habitat.
# for interactive visual exploration. no saved outputs.
# not working yet - bray-curtis should be performed on the actual OTU table, 
# with sequence depth bin as a label. habitat can be a polygon or label to (hopefully) cluster samples.
# ideally samples WILL cluster in similarity by habitat, and not by sequence depth.

source("paths.r")

# load packages
library(vegan)
library(colorRampPalette)
library(RColorBrewer)


#----------------------------- PRIORS - BAHRAM -----------------------------------------#

### data prep ###

# read in descriptors for samples & sites 
map<-readRDS("/fs/data3/caverill/NEFI_data/16S/pecan_gen/bahram_prior_gen.rds")
map <- readRDS("16S/data_construction/data_QC/metadata_all_lats.rds")
# read in OTU file
bahram_otu <- readRDS(bahram_dada2_SV_table.path)
#bahram_otu <- readRDS(bahram_dada2_SV_table_rare.path)
map <- as.data.frame(map)
map <- map[as.character(map$Run) %in% rownames(bahram_otu),]
# subset by sites in dataset
bahram_otu <- bahram_otu[rownames(bahram_otu) %in% map$Run,]

# read in sequence depth files
bahram_track <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/bahram_processed/track.rds")
bahram_track <- as.data.frame(bahram_track)
bahram_track$Run <- rownames(bahram_track)
# subset by sites in dataset
bahram_track <- bahram_track[rownames(bahram_track) %in% map$Run,]

# get covariates for plot
Seq_bin <- as.numeric(cut(bahram_track$seqs_nonchim, 12))
treat <- map$Biome
levels(treat)[levels(treat)=="-"] <- "other"
PH <- cut(map$pH, breaks=8)
CN <- cut(map$cn, breaks=8)
MAP <- cut(map$map, breaks=8)
MAT <- cut(map$mat, breaks=8)
lat <- cut(map$Lat, breaks=8)
### visualization ###

# BRAY-CURTIS #

# set up matrices #
dist.matrix<-as.matrix(bahram_otu)
init_NMDS=metaMDS(dist.matrix,k=2,trymax=1000,distance = "bray")
all.site.coordinates<-data.frame(scores(init_NMDS))

# view stressplot
stressplot(init_NMDS)

# view sample names on plot - not that useful
ordiplot(init_NMDS,type="n")
orditorp(init_NMDS,display="sites",cex=.6,air=0.01)

# plot samples by seq depth
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - sequencing depth")
points(all.site.coordinates$NMDS1,all.site.coordinates$NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(Seq_bin))],cex=1.5)
#ordiellipse(all.site.coordinates, Seq_bin, kind = "se", draw = "lines", conf = 0.95)
dev.copy(png,'16S/data_construction/data_QC/nmds_seqdepth_bahram.png')
dev.off()
#ordispider(init_NMDS, Seq_bin, col = "purple", label= TRUE)

# plot data by biome 
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - habitat");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(7, "Set2")[unclass(factor(treat))],cex=1.5); 
#ordiellipse(all.site.coordinates, treat, kind = "se", draw = "lines", conf = 0.95)
dev.copy(png,'16S/data_construction/data_QC/nmds_habitat_bahram.png')
dev.off()
#ordispider(init_NMDS, treat, col = "purple",  label= TRUE)

# plot by pH
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - pH");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(PH))],cex=1.5); 


# plot by C:N
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - CN");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(CN))],cex=1.5); 

# plot by MAP
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - MAP");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(MAP))],cex=1.5); 

# plot by MAT
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - MAT");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(MAT))],cex=1.5); 

# plot by lat
plot(init_NMDS,type="n",display = "sites", main="Bahram - Bray-curtis - Latitude");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(lat))],cex=1.5); 


# JACCARD #

# set up matrices #
dist.matrix<-as.matrix(bahram_otu)
jacc_NMDS=metaMDS(dist.matrix,k=2,trymax=1000,distance = "jaccard")
jacc_all.site.coordinates<-data.frame(scores(jacc_NMDS))

# view stressplot
stressplot(jacc_NMDS)

# view sample names on plot - not that useful
ordiplot(jacc_NMDS,type="n")
orditorp(jacc_NMDS,display="sites",cex=.6,air=0.01)

# plot samples by seq depth
plot(jacc_NMDS,type="n",display = "sites")
points(jacc_all.site.coordinates$NMDS1,jacc_all.site.coordinates$NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(Seq_bin))],cex=1.5)
ordiellipse(jacc_all.site.coordinates, Seq_bin, kind = "se", draw = "lines", conf = 0.95)
ordispider(jacc_NMDS, Seq_bin, col = "purple", label= TRUE)

# plot data by treatment 
plot(jacc_NMDS,type="n",display = "sites");
points(jacc_all.site.coordinates$NMDS1, jacc_all.site.coordinates$NMDS2,pch=c(21),bg=brewer.pal(7, "Set2")[unclass(factor(treat))],cex=1.5); 
#ordiellipse(jacc_all.site.coordinates, treat, kind = "se", draw = "lines", conf = 0.95)
ordispider(jacc_NMDS, treat, col = "purple",  label= TRUE)



#----------------------------- NEON -----------------------------------------#


### data prep ###

# read in descriptors for samples & sites 
map<-readRDS("/fs/data3/caverill/NEFI_data/16S/pecan_gen/NEON_data_aggregation/obs.table_16S.rds")
hm <- read.csv("/fs/data3/caverill/NEFI_data/16S/pecan_gen/NEON_covariates/field-sites.csv")
map <- merge(map,hm[,c("Site.ID","Dominant.NLCD.Classes")], by.x="siteID", by.y="Site.ID")

# read in OTU file
NEON_otu <- readRDS(NEON_dada2_SV_table.path)
# subset by sites in dataset
NEON_otu <- NEON_otu[rownames(NEON_otu) %in% map$deprecatedVialID,]

# read in sequence depth files
neon_track <- readRDS("/fs/data3/caverill/NEFI_data/16S/scc_gen/NEON_processed/track.rds")
neon_track <- as.data.frame(neon_track)
# subset by sites in dataset
neon_track <- neon_track[rownames(neon_track) %in% map$deprecatedVialID,]

# remove sites with less than 1k reads
NEON_otu <- NEON_otu[-which(rowSums(NEON_otu)<1000),]
neon_track <- neon_track[rownames(neon_track) %in% rownames(NEON_otu),]
map <- map[map$deprecatedVialID %in% rownames(NEON_otu),]

# WHY TF ARE THERE DUPLICATES. 4 samples were sequenced twice? 3 of those, at both Argonne and Batelle?
dup <- map$deprecatedVialID[duplicated(map$deprecatedVialID)]
repeats <- map[map$deprecatedVialID %in% dup,]
repeats <- repeats[order(repeats$deprecatedVialID),]
# REMOVE THEM. data is the same within this df anyways.
map <- map[-which(duplicated(map$deprecatedVialID)),]


# get covariates for plot
Seq_bin <- as.numeric(cut(neon_track$seqs_nonchim, 12))
treat <- map$Dominant.NLCD.Classes


### visualization ###

# BRAY-CURTIS #

# set up matrices #
dist.matrix<-as.matrix(NEON_otu)
init_NMDS=metaMDS(dist.matrix,k=2,trymax=1000,distance = "bray")
all.site.coordinates<-data.frame(scores(init_NMDS))

# view stressplot
stressplot(init_NMDS)

# view sample names on plot - not that useful
ordiplot(init_NMDS,type="n")
orditorp(init_NMDS,display="sites",cex=.6,air=0.01)

# plot samples by seq depth
plot(init_NMDS,type="n",display = "sites", main="NEON - Bray-curtis - sequencing depth")
points(all.site.coordinates$NMDS1,all.site.coordinates$NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(Seq_bin))],cex=1.5)
dev.copy(png,'16S/data_construction/data_QC/nmds_seqdepth_NEON.png')
dev.off()
#ordiellipse(all.site.coordinates, Seq_bin, kind = "se", draw = "lines", conf = 0.95)
ordispider(init_NMDS, Seq_bin, col = "purple", label= TRUE)

# plot data by treatment 
plot(init_NMDS,type="n",display = "sites", main="NEON - Bray-curtis - habitat");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=brewer.pal(7, "Set2")[unclass(factor(treat))],cex=1.5); 
dev.copy(png,'16S/data_construction/data_QC/nmds_habitat_NEON.png')
dev.off()
#ordiellipse(all.site.coordinates, treat, kind = "se", draw = "lines", conf = 0.95)
ordispider(init_NMDS, treat, col = "purple",  label= F)

# JACCARD #

# set up matrices #
dist.matrix<-as.matrix(NEON_otu)
jacc_NMDS=metaMDS(dist.matrix,k=2,trymax=1000,distance = "jaccard")
jacc_all.site.coordinates<-data.frame(scores(jacc_NMDS))

# view stressplot
stressplot(jacc_NMDS)

# view sample names on plot - not that useful
ordiplot(jacc_NMDS,type="n")
orditorp(jacc_NMDS,display="sites",cex=.6,air=0.01)

# plot samples by seq depth
plot(jacc_NMDS,type="n",display = "sites")
points(jacc_all.site.coordinates$NMDS1,jacc_all.site.coordinates$NMDS2,pch=c(21),bg=brewer.pal(9, "RdPu")[unclass(factor(Seq_bin))],cex=1.5)
ordiellipse(jacc_all.site.coordinates, Seq_bin, kind = "se", draw = "lines", conf = 0.95)
ordispider(jacc_NMDS, Seq_bin, col = "purple", label= TRUE)

# plot data by treatment 
plot(jacc_NMDS,type="n",display = "sites");
points(jacc_all.site.coordinates$NMDS1, jacc_all.site.coordinates$NMDS2,pch=c(21),bg=brewer.pal(7, "Set2")[unclass(factor(treat))],cex=1.5); 
#ordiellipse(jacc_all.site.coordinates, treat, kind = "se", draw = "lines", conf = 0.95)
ordispider(jacc_NMDS, treat, col = "purple",  label= TRUE)






















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

all.MDS<-metaMDS(as.dist(log(dist.matrix +1)),k=3, trymax = 150,zerodist="add")

#set up plot
all.site.coordinates<-data.frame(scores(all.MDS))

#plot samples with Sample ID
ordiplot(all.MDS,type="text")
ordihull(all.MDS,groups=treat,draw="polygon",col="grey90",label=F)

#plot samples by seq depth
plot(all.MDS,type="n",display = "sites")
points(all.site.coordinates$NMDS1,all.site.coordinates$NMDS2,pch=c(21),bg=c("aliceblue","steelblue1","steelblue4","darkblue")[unclass(factor(Seq_bin))],cex=2.5)
ordiellipse(all.site.coordinates, Seq_bin, kind = "se", draw = "lines", conf = 0.95)
ordispider(all.MDS, Seq_bin, col = "blue", label= TRUE)


#plot data by treatment (color)
plot(all.MDS,type="n",display = "sites");
points(all.site.coordinates $NMDS1, all.site.coordinates $NMDS2,pch=c(21),bg=c("aliceblue","indianred1","turquoise3")[unclass(factor(treat))],cex=2.5); 
ordiellipse(all.site.coordinates, treat, kind = "se", draw = "lines", conf = 0.95)
